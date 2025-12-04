/***********************************************************************
MDVisualizer - Vrui application class to show time-varying MD trace
data.
Copyright (c) 2005-2025 Oliver Kreylos

This file is part of the MD Visualizer (MDVisualizer).

The MD Visualizer is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

The MD Visualizer is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with the MD Visualizer; if not, write to the Free Software Foundation,
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
***********************************************************************/

#include "MDVisualizer.h"

#include <string>
#include <stdexcept>
#include <Misc/StdError.h>
#include <Misc/FileNameExtensions.h>
#include <Misc/MessageLogger.h>
#include <Cluster/MulticastPipe.h>
#include <Math/Math.h>
#include <Geometry/OrthogonalTransformation.h>
#include <GL/gl.h>
#include <GL/GLVertexTemplates.h>
#include <GL/GLFog.h>
#include <GL/GLContextData.h>
#include <GL/GLExtensionManager.h>
#include <GL/Extensions/GLARBVertexBufferObject.h>
#include <GL/GLModels.h>
#include <GLMotif/StyleSheet.h>
#include <GLMotif/PopupMenu.h>
#include <GLMotif/PopupWindow.h>
#include <GLMotif/Margin.h>
#include <GLMotif/RowColumn.h>
#include <GLMotif/Blind.h>
#include <GLMotif/Label.h>
#include <GLMotif/TextField.h>
#include <GLMotif/Button.h>
#include <Vrui/Vrui.h>
#include <Collaboration2/Client.h>

#include "Color.h"
#include "Trace.h"
#include "CharmmTrace.h"
#include "LammpsTrace.h"
#include "EspressoTrace.h"

// DEBUGGING
// #include <Realtime/Time.h>

/*********************************************
Declaration of class MDVisualizer::BondFinder:
*********************************************/

class MDVisualizer::BondFinder
	{
	/* Elements: */
	private:
	const AtomPoint& searchAtom; // Pointer to the atom for which to find bonds
	MD::Scalar maxBondLength; // Maximum bond length involving the search atom
	BondList& bonds; // List to which to add found bonds
	
	/* Constructors and destructors: */
	public:
	BondFinder(const AtomPoint& sSearchAtom,MD::Scalar sMaxHalfBondLength,BondList& sBonds)
		:searchAtom(sSearchAtom),
		 maxBondLength(searchAtom.value.second+sMaxHalfBondLength),
		 bonds(sBonds)
		{
		}
	
	/* Methods: */
	const MD::Point& getQueryPosition(void) const
		{
		return searchAtom;
		}
	bool operator()(const AtomPoint& node,int splitDimension)
		{
		/* Check if the node atom can bond with the search atom: */
		if(searchAtom.value.first<node.value.first) // Check indices to get each bond only once
			{
			MD::Scalar bondLength2=Math::sqr(searchAtom.value.second+node.value.second);
			if(Geometry::sqrDist(searchAtom,node)<=bondLength2)
				{
				/* Add the bond to the list: */
				bonds.push_back(Bond(searchAtom.value.first,node.value.first));
				}
			}
		
		/* Check if the other half of the split plane can be culled: */
		return Math::abs(node[splitDimension]-searchAtom[splitDimension])<=maxBondLength;
		}
	};

/****************************************************
Declaration of struct MDVisualizer::BondFinderThread:
****************************************************/

struct MDVisualizer::BondFinderThread
	{
	/* Elements: */
	public:
	bool workRequest; // Flag if a work request is pending
	const AtomTree* atomTree; // Pointer to kd-tree containing atom records
	const AtomPoint* atomsBegin; // Pointer to beginning of sub-array of atom records in a kd-tree
	const AtomPoint* atomsEnd; // Pointer to end of sub-array of atom records in a kd-tree
	Threads::Thread thread; // The thread executing the bond finding algorithm
	BondList bonds; // List of bonds found in the sub-array of atom records
	};

/***************************************
Methods of class MDVisualizer::DataItem:
***************************************/

MDVisualizer::DataItem::DataItem(void)
	:bufferVersion(0),
	 influenceSphereDisplayListId(0)
	{
	/* Initialize the vertex buffer object extension: */
	GLARBVertexBufferObject::initExtension();
	
	/* Create buffer objects: */
	glGenBuffersARB(1,&atomVertexBufferObjectId);
	glGenBuffersARB(1,&bondIndexBufferObjectId);
	
	/* Create display lists: */
	influenceSphereDisplayListId=glGenLists(1);
	}

MDVisualizer::DataItem::~DataItem(void)
	{
	/* Delete vertex buffer objects: */
	glDeleteBuffersARB(1,&atomVertexBufferObjectId);
	glDeleteBuffersARB(1,&bondIndexBufferObjectId);
	
	/* Delete display lists: */
	glDeleteLists(influenceSphereDisplayListId,1);
	}

/*****************************
Methods of class MDVisualizer:
*****************************/

void* MDVisualizer::bondFinderThreadMethod(MDVisualizer::BondFinderThread* bondFinderThread)
	{
	/* Synchronize on the barrier to make sure all threads are started: */
	bondFindersDone.synchronize();
	
	/* Wait for bond finding requests until shut down: */
	while(true)
		{
		/* Wait for the next bond finding request: */
		const AtomTree* atomTree;
		const AtomPoint* atomsBegin;
		const AtomPoint* atomsEnd;
		{
		Threads::MutexCond::Lock bondFinderRequestLock(bondFinderRequestCond);
		while(bondFinderThreadsRunning&&!bondFinderThread->workRequest)
			bondFinderRequestCond.wait(bondFinderRequestLock);
		atomTree=bondFinderThread->atomTree;
		atomsBegin=bondFinderThread->atomsBegin;
		atomsEnd=bondFinderThread->atomsEnd;
		bondFinderThread->workRequest=false;
		}
		
		/* Bail out if shutting down: */
		if(!bondFinderThreadsRunning)
			break;
		
		/* Find all bonds between close atoms in the assigned sub-array: */
		{
		bondFinderThread->bonds.clear();
		for(const AtomPoint* apPtr=atomsBegin;apPtr!=atomsEnd;++apPtr)
			{
			/* Create a bond finder object for this atom and traverse the atom tree: */
			BondFinder bf(*apPtr,maxHalfBondLength,bondFinderThread->bonds);
			bondFinderThread->atomTree->traverseTreeDirected(bf);
			}
		}
		
		/* Signal completion to the loader thread: */
		bondFindersDone.synchronize();
		}
	
	return 0;
	}

void* MDVisualizer::loaderThreadMethod(void)
	{
	/* Wait for time step load requests until shut down: */
	while(true)
		{
		/* Wait for the next load request: */
		size_t timeStepIndex=0;
		{
		Threads::MutexCond::Lock loadRequestLock(loadRequestCond);
		while(loaderThreadRunning&&!loadRequest)
			loadRequestCond.wait(loadRequestLock);
		timeStepIndex=loadTimeStepIndex;
		loadRequest=false;
		}
		
		/* Bail out if shutting down: */
		if(!loaderThreadRunning)
			break;
		
		// DEBUGGING
		// Realtime::TimePointMonotonic now;
		
		/* Start a new slot in the time step triple buffer: */
		TimeStep& timeStep=timeSteps.startNewValue();
		
		/* Load the requested time step: */
		timeStep.index=timeStepIndex;
		
		/* Retrieve the number of atoms in the time step: */
		size_t newNumAtoms=trace->getNumAtoms(timeStepIndex);
		
		/* Check the size of the atoms array: */
		if(timeStep.numAtoms!=newNumAtoms)
			{
			delete[] timeStep.atoms;
			timeStep.numAtoms=newNumAtoms;
			timeStep.atoms=new MD::Atom[timeStep.numAtoms];
			}
		
		try
			{
			/* Load the time step's atoms: */
			trace->loadAtoms(timeStepIndex,timeStep.atoms);
			}
		catch(const std::runtime_error& err)
			{
			Misc::formattedUserError("MDVisualizer::loaderThread: Unable to load time step %u due to exception %s",(unsigned int)(timeStepIndex),err.what());
			continue;
			}
		
		/* Create the atom kd-tree: */
		AtomPoint* atomPoints=timeStep.atomTree.createTree(timeStep.numAtoms);
		AtomPoint* apPtr=atomPoints;
		MD::Atom* aPtr=timeStep.atoms;
		for(unsigned int i=0;i<timeStep.numAtoms;++i,++apPtr,++aPtr)
			{
			*apPtr=aPtr->getPosition();
			apPtr->value.first=i;
			apPtr->value.second=aPtr->getHalfBondLength(*atomTypeMap);
			}
		timeStep.atomTree.releasePoints();
		
		/* Delegate work to the bond finder threads: */
		{
		Threads::MutexCond::Lock bondFinderRequestLock(bondFinderRequestCond);
		for(unsigned int i=0;i<numBondFinderThreads;++i)
			{
			bondFinderThreads[i].atomTree=&timeStep.atomTree;
			bondFinderThreads[i].atomsBegin=atomPoints+(((unsigned int)(timeStep.numAtoms)*i)/numBondFinderThreads);
			bondFinderThreads[i].atomsEnd=atomPoints+(((unsigned int)(timeStep.numAtoms)*(i+1))/numBondFinderThreads);
			bondFinderThreads[i].workRequest=true;
			}
		bondFinderRequestCond.broadcast();
		}
		
		/* Wait until all bond finder threads have finished: */
		bondFindersDone.synchronize();
		
		/* Assemble the results from all bond finder threads: */
		timeStep.bonds.clear();
		for(unsigned int i=0;i<numBondFinderThreads;++i)
			timeStep.bonds.insert(timeStep.bonds.end(),bondFinderThreads[i].bonds.begin(),bondFinderThreads[i].bonds.end());
		
		// DEBUGGING
		// double time(now.setAndDiff());
		// Misc::formattedLogNote("MDVisualizer::loaderThread: Loaded %u atoms and created %u bonds in %f ms",timeStep.numAtoms,timeStep.bonds.size(),time*1000.0);
		
		/* Post the new time step: */
		timeSteps.postNewValue();
		Vrui::requestUpdate();
		}
	
	return 0;
	}

void MDVisualizer::loadTimeStep(size_t newTimeStepIndex)
	{
	Threads::MutexCond::Lock loadRequestLock(loadRequestCond);
	loadRequest=true;
	loadTimeStepIndex=newTimeStepIndex;
	loadRequestCond.signal();
	}

void MDVisualizer::loadTimeStepIndexUpdatedCallback(KoinoniaClient* client,KoinoniaProtocol::ObjectID id,void* object,void* userData)
	{
	MDVisualizer* thisPtr=static_cast<MDVisualizer*>(userData);
	
	/* Request loading the new time step: */
	thisPtr->loadTimeStep(thisPtr->loadTimeStepIndex);
	}

void MDVisualizer::renderingParametersUpdatedCallback(KoinoniaClient* client,KoinoniaProtocol::ObjectID id,void* object,void* userData)
	{
	MDVisualizer* thisPtr=static_cast<MDVisualizer*>(userData);
	
	/* Update the UI: */
	thisPtr->renderDialog->updateVariables();
	
	#if USE_IMPOSTORS
	
	/* Update the atom and bond renderers: */
	thisPtr->sphereRenderer.setFixedRadius(thisPtr->renderingParameters.atomRadius);
	thisPtr->cylinderRenderer.setFixedRadius(thisPtr->renderingParameters.bondRadius);
	
	#endif
	}

void MDVisualizer::clearSelectionCallback(Misc::CallbackData* cbData)
	{
	/* Clear the set of selected atoms: */
	selectedAtomsSet.clear();
	selectedAtoms.clear();
	}

void MDVisualizer::showAnimationDialogCallback(Misc::CallbackData* cbData)
	{
	/* Pop up the animation dialog: */
	Vrui::popupPrimaryWidget(animationDialog);
	}

void MDVisualizer::showRenderDialogCallback(Misc::CallbackData* cbData)
	{
	/* Pop up the render dialog: */
	Vrui::popupPrimaryWidget(renderDialog);
	}

GLMotif::PopupMenu* MDVisualizer::createMainMenu(void)
	{
	GLMotif::PopupMenu* mainMenu=new GLMotif::PopupMenu("MainMenuPopup",Vrui::getWidgetManager());
	mainMenu->setTitle("MD Trace Visualizer");
	
	mainMenu->addEntry("Clear Selection")->getSelectCallbacks().add(this,&MDVisualizer::clearSelectionCallback);
	
	GLMotif::Button* showAnimationDialogButton=new GLMotif::Button("ShowAnimationDialogButton",mainMenu,"Show Animation Dialog");
	showAnimationDialogButton->getSelectCallbacks().add(this,&MDVisualizer::showAnimationDialogCallback);
	
	GLMotif::Button* showRenderDialogButton=new GLMotif::Button("ShowRenderDialogButton",mainMenu,"Show Render Dialog");
	showRenderDialogButton->getSelectCallbacks().add(this,&MDVisualizer::showRenderDialogCallback);
	
	mainMenu->manageMenu();
	return mainMenu;
	}

void MDVisualizer::timeStepIndexCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	size_t newTimeStepIndex=size_t(Math::floor(cbData->value+0.5));
	
	/* Request loading the new time step: */
	if(newTimeStepIndex!=loadTimeStepIndex)
		loadTimeStep(newTimeStepIndex);
	
	/* Share the new time step index with the server: */
	if(koinonia!=0)
		koinonia->replaceSharedObject(loadTimeStepIndexId);
	}

void MDVisualizer::animationSpeedCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData)
	{
	animationSpeed=cbData->value;
	}

void MDVisualizer::animateCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData)
	{
	if(cbData->set)
		{
		/* Start animating: */
		animate=true;
		nextTimeStepTime=Vrui::getApplicationTime()+1.0/animationSpeed;
		
		/* Request loading the next time step: */
		size_t newTimeStepIndex=loadTimeStepIndex+1;
		if(newTimeStepIndex>=trace->getNumTimeSteps())
			newTimeStepIndex=0;
		loadTimeStep(newTimeStepIndex);
		
		/* Share the new time step index with the server: */
		if(koinonia!=0)
			koinonia->replaceSharedObject(loadTimeStepIndexId);
		
		/* Update the UI: */
		timeStepIndexSlider->setValue(loadTimeStepIndexId);
		}
	else
		{
		/* Stop animating: */
		animate=false;
		nextTimeStepTime=0.0;
		}
	}

GLMotif::PopupWindow* MDVisualizer::createAnimationDialog(void)
	{
	const GLMotif::StyleSheet& ss=*Vrui::getUiStyleSheet();
	
	GLMotif::PopupWindow* animationDialogPopup=new GLMotif::PopupWindow("AnimationDialogPopup",Vrui::getWidgetManager(),"Animation Dialog");
	
	GLMotif::RowColumn* animation=new GLMotif::RowColumn("Animation",animationDialogPopup,false);
	animation->setNumMinorWidgets(2);
	
	new GLMotif::Label("TraceFileIndexLabel",animation,"Trace File Index");
	
	timeStepIndexSlider=new GLMotif::TextFieldSlider("TimeStepIndexSlider",animation,7,ss.fontHeight*20.0f);
	timeStepIndexSlider->setSliderMapping(GLMotif::TextFieldSlider::LINEAR);
	timeStepIndexSlider->setValueType(GLMotif::TextFieldSlider::UINT);
	timeStepIndexSlider->setValueRange(0,trace->getNumTimeSteps()-1,1.0);
	timeStepIndexSlider->setValue(loadTimeStepIndex);
	timeStepIndexSlider->getValueChangedCallbacks().add(this,&MDVisualizer::timeStepIndexCallback);
	
	new GLMotif::Label("AnimationSpeedLabel",animation,"Animation Speed");
	
	GLMotif::TextFieldSlider* animationSpeedSlider=new GLMotif::TextFieldSlider("AnimationSpeedSlider",animation,7,ss.fontHeight*20.0f);
	animationSpeedSlider->setSliderMapping(GLMotif::TextFieldSlider::LINEAR);
	animationSpeedSlider->setValueType(GLMotif::TextFieldSlider::FLOAT);
	animationSpeedSlider->setValueRange(0.1,100.0,0.1);
	animationSpeedSlider->setValue(animationSpeed);
	animationSpeedSlider->getValueChangedCallbacks().add(this,&MDVisualizer::animationSpeedCallback);
	
	new GLMotif::Blind("Blind1",animation);
	
	GLMotif::Margin* toggleMargin=new GLMotif::Margin("ToggleMargin",animation,false);
	toggleMargin->setAlignment(GLMotif::Alignment::LEFT);
	
	animateToggle=new GLMotif::ToggleButton("AnimateToggle",toggleMargin,"Animate");
	animateToggle->setBorderType(GLMotif::Widget::PLAIN);
	animateToggle->setBorderWidth(0.0f);
	animateToggle->setToggle(animate);
	animateToggle->getValueChangedCallbacks().add(this,&MDVisualizer::animateCallback);
	
	toggleMargin->manageChild();
	
	animation->manageChild();
	
	return animationDialogPopup;
	}

void MDVisualizer::renderingParametersChangedCallback(Misc::CallbackData* cbData)
	{
	/* Share the new rendering parameters with the server: */
	if(koinonia!=0)
		koinonia->replaceSharedObject(renderingParametersId);
	
	#if USE_IMPOSTORS
	
	/* Update the atom and bond renderers: */
	sphereRenderer.setFixedRadius(renderingParameters.atomRadius);
	cylinderRenderer.setFixedRadius(renderingParameters.bondRadius);
	
	#endif
	}

GLMotif::PopupWindow* MDVisualizer::createRenderDialog(void)
	{
	const GLMotif::StyleSheet& ss=*Vrui::getUiStyleSheet();
	
	GLMotif::PopupWindow* renderDialogPopup=new GLMotif::PopupWindow("RenderDialogPopup",Vrui::getWidgetManager(),"Render Dialog");
	
	GLMotif::RowColumn* render=new GLMotif::RowColumn("Render",renderDialogPopup,false);
	render->setOrientation(GLMotif::RowColumn::VERTICAL);
	render->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	render->setNumMinorWidgets(2);
	
	new GLMotif::Blind("Blind1",render);
	
	GLMotif::Margin* toggleMargin=new GLMotif::Margin("ToggleMargin",render,false);
	toggleMargin->setAlignment(GLMotif::Alignment::LEFT);
	
	GLMotif::RowColumn* toggleBox=new GLMotif::RowColumn("ToggleBox",toggleMargin,false);
	toggleBox->setOrientation(GLMotif::RowColumn::HORIZONTAL);
	toggleBox->setPacking(GLMotif::RowColumn::PACK_TIGHT);
	toggleBox->setNumMinorWidgets(1);
	
	GLMotif::ToggleButton* drawAtomsToggle=new GLMotif::ToggleButton("DrawAtomsToggle",toggleBox,"Draw Atoms");
	drawAtomsToggle->setBorderType(GLMotif::Widget::PLAIN);
	drawAtomsToggle->setBorderWidth(0.0f);
	drawAtomsToggle->track(renderingParameters.drawAtoms);
	drawAtomsToggle->getValueChangedCallbacks().add(this,&MDVisualizer::renderingParametersChangedCallback);
	
	GLMotif::ToggleButton* drawBondsToggle=new GLMotif::ToggleButton("DrawBondsToggle",toggleBox,"Draw Bonds");
	drawBondsToggle->setBorderType(GLMotif::Widget::PLAIN);
	drawBondsToggle->setBorderWidth(0.0f);
	drawBondsToggle->track(renderingParameters.drawBonds);
	drawBondsToggle->getValueChangedCallbacks().add(this,&MDVisualizer::renderingParametersChangedCallback);
	
	toggleBox->manageChild();
	
	toggleMargin->manageChild();
	
	new GLMotif::Label("AtomRadiusLabel",render,"Atom Radius");
	
	GLMotif::TextFieldSlider* atomRadiusSlider=new GLMotif::TextFieldSlider("AtomRadiusSlider",render,6,ss.fontHeight*10.0f);
	atomRadiusSlider->setSliderMapping(GLMotif::TextFieldSlider::LINEAR);
	atomRadiusSlider->setValueType(GLMotif::TextFieldSlider::FLOAT);
	atomRadiusSlider->setValueRange(0.01,1.0,0.01);
	atomRadiusSlider->track(renderingParameters.atomRadius);
	atomRadiusSlider->getValueChangedCallbacks().add(this,&MDVisualizer::renderingParametersChangedCallback);
	
	new GLMotif::Label("BondRadiusLabel",render,"Bond Radius");
	
	GLMotif::TextFieldSlider* bondRadiusSlider=new GLMotif::TextFieldSlider("BondRadiusSlider",render,6,ss.fontHeight*10.0f);
	bondRadiusSlider->setSliderMapping(GLMotif::TextFieldSlider::LINEAR);
	bondRadiusSlider->setValueType(GLMotif::TextFieldSlider::FLOAT);
	bondRadiusSlider->setValueRange(0.01,1.0,0.01);
	bondRadiusSlider->track(renderingParameters.bondRadius);
	bondRadiusSlider->getValueChangedCallbacks().add(this,&MDVisualizer::renderingParametersChangedCallback);
	
	render->manageChild();
	
	return renderDialogPopup;
	}

MDVisualizer::MDVisualizer(int& argc,char**& argv)
	:Vrui::Application(argc,argv),
	 trace(0),
	 koinonia(0),
	 loadRequest(false),loadTimeStepIndex(0),
	 timeStepPipe(Vrui::openPipe()),
	 loaderThreadRunning(false),
	 bondFinderThreadsRunning(false),numBondFinderThreads(4),bondFinderThreads(new BondFinderThread[numBondFinderThreads]),
	 bondFindersDone(numBondFinderThreads+1),
	 timeStepVersion(0),haveTimeStep(false),
	 selectedAtomsSet(101),
	 atomMaterial(GLMaterial::Color(0.333f,0.333f,0.333f),GLMaterial::Color(1.0f,1.0f,1.0f),25.0f),
	 animationSpeed(10.0),animate(false),nextTimeStepTime(0.0),
	 mainMenu(0),animationDialog(0),renderDialog(0)
	{
	/* Parse the command line: */
	int traceType=-1;
	const char* traceFileName=0;
	for(int argi=1;argi<argc;++argi)
		{
		if(argv[argi][0]=='-')
			{
			if(strcasecmp(argv[argi]+1,"charmm")==0)
				traceType=0;
			else if(strcasecmp(argv[argi]+1,"lammps")==0)
				traceType=1;
			else if(strcasecmp(argv[argi]+1,"espresso")==0)
				traceType=2;
			else
				Misc::formattedConsoleWarning("MDVisualizer::MDVisualizer: Ignoring command line option %s",argv[argi]);
			}
		else if(traceFileName==0)
			traceFileName=argv[argi];
		else
			Misc::formattedConsoleWarning("MDVisualizer::MDVisualizer: Ignoring command line argument %s",argv[argi]);
		}
	if(traceFileName==0)
		throw std::runtime_error("MDVisualizer::MDVisualizer: No trace file name provided");
	if(traceType==-1)
		{
		if(Misc::hasCaseExtension(traceFileName,".crd"))
			traceType=0;
		else if(Misc::hasCaseExtension(traceFileName,".lammpstrj"))
			traceType=1;
		else
			throw Misc::makeStdErr(__PRETTY_FUNCTION__,"Unrecognized trace file format");
		}
	if(traceType==0)
		{
		/* Load a CHARMM trace file: */
		trace=new MD::CharmmTrace(traceFileName);
		}
	else if(traceType==1)
		{
		/* Construct a meta data file name: */
		std::string metaDataFileName=traceFileName;
		metaDataFileName.append(".meta");
		
		/* Load a LAMMP trace file: */
		trace=new MD::LammpsTrace(metaDataFileName.c_str(),traceFileName);
		}
	else
		{
		/* Load an Espresso trace file: */
		trace=new MD::EspressoTrace(traceFileName);
		}
	
	/* Retrieve the trace file's atom type map, bounding box, and maximum half bond length: */
	atomTypeMap=&trace->getAtomTypeMap();
	boundingBox=trace->calcBoundingBox();
	maxHalfBondLength=trace->calcMaxHalfBondLength();
	
	/* Start the bond finder threads: */
	bondFinderThreadsRunning=true;
	for(unsigned int i=0;i<numBondFinderThreads;++i)
		{
		bondFinderThreads[i].workRequest=false;
		bondFinderThreads[i].atomTree=0;
		bondFinderThreads[i].atomsBegin=0;
		bondFinderThreads[i].atomsEnd=0;
		bondFinderThreads[i].thread.start(this,&MDVisualizer::bondFinderThreadMethod,bondFinderThreads+i);
		}
	
	/* Synchronize on the barrier to make sure all threads are started: */
	bondFindersDone.synchronize();
	
	/* Start the trace loader thread: */
	loaderThreadRunning=true;
	loaderThread.start(this,&MDVisualizer::loaderThreadMethod);
	
	/* Initialize rendering parameters: */
	renderingParameters.drawAtoms=true;
	renderingParameters.atomRadius=0.3f;
	renderingParameters.drawBonds=true;
	renderingParameters.bondRadius=0.2f;
	renderingParameters.selectedAtomRadius=0.5f;
	
	#if USE_IMPOSTORS
	
	/* Initialize the sphere and cylinder renderers: */
	sphereRenderer.setFixedRadius(renderingParameters.atomRadius);
	sphereRenderer.setColorMaterial(true);
	cylinderRenderer.setFixedRadius(renderingParameters.bondRadius);
	cylinderRenderer.setCapped(false);
	cylinderRenderer.setColorMaterial(true);
	cylinderRenderer.setBicolor(true);
	selectedSphereRenderer.setVariableRadius();
	selectedSphereRenderer.setColorMaterial(true);
	
	#endif
	
	/* Create the user interface: */
	mainMenu=createMainMenu();
	Vrui::setMainMenu(mainMenu);
	animationDialog=createAnimationDialog();
	renderDialog=createRenderDialog();
	
	/* Check if there is a collaboration client: */
	Collab::Client* client=Collab::Client::getTheClient();
	if(client!=0)
		{
		/* Request a Koinonia client to share the array of enabled flags: */
		koinonia=static_cast<KoinoniaClient*>(client->requestPluginProtocol("Koinonia"));
		
		/* Share the load time step index and the rendering parameters: */
		DataType dataType;
		
		DataType::TypeID loadTimeStepIndexTypeId=DataType::UInt64;
		
		DataType::StructureElement renderingParametersElems[]=
			{
			{DataType::Bool,offsetof(RenderingParameters,drawAtoms)},
			{DataType::Float32,offsetof(RenderingParameters,atomRadius)},
			{DataType::Bool,offsetof(RenderingParameters,drawBonds)},
			{DataType::Float32,offsetof(RenderingParameters,bondRadius)},
			{DataType::Float32,offsetof(RenderingParameters,selectedAtomRadius)}
			};
		DataType::TypeID renderingParametersTypeId=dataType.createStructure(5,renderingParametersElems,sizeof(RenderingParameters));
		
		loadTimeStepIndexId=koinonia->shareObject("MDVisualizer.loadTimeStepIndex",(1U<<16)+0U,dataType,loadTimeStepIndexTypeId,&loadTimeStepIndex,&MDVisualizer::loadTimeStepIndexUpdatedCallback,this);
		renderingParametersId=koinonia->shareObject("MDVisualizer.renderingParameters",(1U<<16)+0U,dataType,renderingParametersTypeId,&renderingParameters,&MDVisualizer::renderingParametersUpdatedCallback,this);
		}
	else
		{
		/* Immediately load the first time step: */
		loadTimeStep(0);
		}
	}

MDVisualizer::~MDVisualizer(void)
	{
	/* Shut down the trace loader thread: */
	loaderThreadRunning=false;
	loadRequestCond.signal();
	loaderThread.join();
	delete timeStepPipe;
	
	/* Shut down the bond finder threads: */
	bondFinderThreadsRunning=false;
	bondFinderRequestCond.broadcast();
	for(unsigned int i=0;i<numBondFinderThreads;++i)
		bondFinderThreads[i].thread.join();
	delete[] bondFinderThreads;
	
	/* Clean up: */
	delete trace;
	delete mainMenu;
	delete animationDialog;
	delete renderDialog;
	}

void MDVisualizer::frame(void)
	{
	/* Check if a new trace time step has arrived and the animation time has passed: */
	if(Vrui::getApplicationTime()>=nextTimeStepTime&&timeSteps.lockNewValue())
		{
		/* Invalidate render buffers: */
		++timeStepVersion;
		
		/* Check if this was the first time step: */
		if(!haveTimeStep)
			{
			/* Initialize the navigation transformation: */
			haveTimeStep=true;
			resetNavigation();
			}
		
		/* Check if animation is running: */
		if(animate)
			{
			/* Load the next time step unless at the end of the trace file: */
			size_t newTimeStepIndex=loadTimeStepIndex+1;
			if(newTimeStepIndex<trace->getNumTimeSteps())
				{
				/* Request loading the next time step: */
				loadTimeStep(newTimeStepIndex);
				
				/* Share the new time step index with the server: */
				if(koinonia!=0)
					koinonia->replaceSharedObject(loadTimeStepIndexId);
				
				/* Update the UI: */
				timeStepIndexSlider->setValue(newTimeStepIndex);
				
				/* Advance the animation timer: */
				nextTimeStepTime+=1.0/animationSpeed;
				}
			else
				{
				/* Stop animation: */
				animate=false;
				nextTimeStepTime=0.0;
				animateToggle->setToggle(false);
				}
			}
		}
	
	if(animate)
		Vrui::scheduleUpdate(nextTimeStepTime);
	}

namespace {

/****************
Helper functions:
****************/

void glColor(const MD::Color& color)
	{
	glColor3fv(color.getRgb());
	}

}

void MDVisualizer::display(GLContextData& contextData) const
	{
	/* Bail out if there is no current time step: */
	if(!haveTimeStep)
		return;
	
	/* Retrieve the context data item: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	/* Render the locked time step buffer: */
	const TimeStep& timeStep=timeSteps.getLockedValue();
	
	/* Set up OpenGL state: */
	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_ENABLE_BIT|GL_FOG_BIT|GL_POINT_BIT|GL_LINE_BIT);
	#if !USE_IMPOSTORS
	glDisable(GL_LIGHTING);
	glPointSize(3.0f);
	glLineWidth(1.0f);
	#endif
	
	/* Enable simplistic depth cueing: */
	glEnable(GL_FOG);
	glFog(GLFog(GLFogEnums::LINEAR,Vrui::getFrontplaneDist(),Vrui::getFrontplaneDist()+Vrui::getDisplaySize()*Vrui::Scalar(6),Vrui::getBackgroundColor()));
	
	GLVertexArrayParts::enable(DataItem::Vertex::getPartsMask());
	glBindBufferARB(GL_ARRAY_BUFFER_ARB,dataItem->atomVertexBufferObjectId);
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,dataItem->bondIndexBufferObjectId);
	
	/* Check if the buffer objects have been invalidated: */
	if(dataItem->bufferVersion!=timeStepVersion)
		{
		/* Upload current trace data to buffer objects: */
		glBufferDataARB(GL_ARRAY_BUFFER_ARB,timeStep.numAtoms*sizeof(DataItem::Vertex),0,GL_DYNAMIC_DRAW_ARB);
		DataItem::Vertex* vertexPtr=static_cast<DataItem::Vertex*>(glMapBufferARB(GL_ARRAY_BUFFER_ARB,GL_WRITE_ONLY_ARB));
		const MD::Atom* atomsEnd=timeStep.atoms+timeStep.numAtoms;
		for(const MD::Atom* atomPtr=timeStep.atoms;atomPtr!=atomsEnd;++atomPtr,++vertexPtr)
			{
			const MD::Atom::Type& atomType=atomPtr->getType(*atomTypeMap);
			for(int i=0;i<3;++i)
				vertexPtr->color[i]=GLubyte(Math::floor(atomType.color[i]*255.0f+0.5f));
			for(int i=0;i<3;++i)
				vertexPtr->position[i]=GLfloat(atomPtr->getPosition()[i]);
			}
		glUnmapBufferARB(GL_ARRAY_BUFFER_ARB);
		
		glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB,timeStep.bonds.size()*2*sizeof(GLuint),0,GL_DYNAMIC_DRAW_ARB);
		GLuint* indexPtr=static_cast<GLuint*>(glMapBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,GL_WRITE_ONLY_ARB));
		for(BondList::const_iterator bIt=timeStep.bonds.begin();bIt!=timeStep.bonds.end();++bIt,indexPtr+=2)
			{
			indexPtr[0]=bIt->first;
			indexPtr[1]=bIt->second;
			}
		glUnmapBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB);
		
		/* Validate the buffer objects: */
		dataItem->bufferVersion=timeStepVersion;
		}
	
	/* Draw all atoms and/or bonds: */
	glMaterial(GLMaterialEnums::FRONT,atomMaterial);
	glVertexPointer(static_cast<DataItem::Vertex*>(0));
	
	if(renderingParameters.drawAtoms)
		{
		/* Draw all atoms: */
		#if USE_IMPOSTORS
		sphereRenderer.enable(GLfloat(Vrui::getNavigationTransformation().getScaling()),contextData);
		glDrawArrays(GL_POINTS,0,timeStep.numAtoms);
		sphereRenderer.disable(contextData);
		#else
		glDrawArrays(GL_POINTS,0,timeStep.numAtoms);
		#endif
		}
	
	if(renderingParameters.drawBonds)
		{
		/* Draw all bonds: */
		#if USE_IMPOSTORS
		cylinderRenderer.enable(GLfloat(Vrui::getNavigationTransformation().getScaling()),contextData);
		glDrawElements(GL_LINES,timeStep.bonds.size()*2,GL_UNSIGNED_INT,static_cast<GLuint*>(0));
		cylinderRenderer.disable(contextData);
		#else
		glDrawElements(GL_LINES,timeStep.bonds.size()*2,GL_UNSIGNED_INT,static_cast<GLuint*>(0));
		#endif
		}
	
	/* Disable vertex arrays and protect the buffers: */
	glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,0);
	GLVertexArrayParts::disable(DataItem::Vertex::getPartsMask());
	
	/* Draw all selected atoms: */
	#if USE_IMPOSTORS
	selectedSphereRenderer.enable(GLfloat(Vrui::getNavigationTransformation().getScaling()),contextData);
	#endif
	glBegin(GL_POINTS);
	for(std::vector<unsigned int>::const_iterator saIt=selectedAtoms.begin();saIt!=selectedAtoms.end();++saIt)
		{
		/* Draw the atom's sphere: */
		const MD::Atom& atom=timeStep.atoms[*saIt];
		const MD::Atom::Type& atomType=atom.getType(*atomTypeMap);
		glColor(atomType.color);
		glVertex(atom.getPosition()[0],atom.getPosition()[1],atom.getPosition()[2],atomType.radius);
		}
	glEnd();
	#if USE_IMPOSTORS
	selectedSphereRenderer.disable(contextData);
	#endif
	
	/* Reset OpenGL state: */
	glPopAttrib();
	}

void MDVisualizer::resetNavigation(void)
	{
	/* Bail out if there is no current time step: */
	if(!haveTimeStep)
		return;
	
	/* Calculate center and radius of a sphere containing all atoms in the current trace: */
	MD::Point center=Geometry::mid(boundingBox.min,boundingBox.max);
	MD::Scalar size=Geometry::dist(boundingBox.min,boundingBox.max);
	
	/* Set the navigation transformation: */
	Vrui::setNavigationTransformation(center,size,Vrui::Vector(0,0,1));
	}

void MDVisualizer::initContext(GLContextData& contextData) const
	{
	/* Create a context data item: */
	DataItem* dataItem=new DataItem;
	contextData.addDataItem(this,dataItem);
	
	/* Create the influence sphere display list: */
	glNewList(dataItem->influenceSphereDisplayListId,GL_COMPILE);
	glDisable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glDepthMask(GL_FALSE);
	glLineWidth(1.0f);
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	glColor4f(1.0f,1.0f,0.0f,0.67f);
	glDrawSphereIcosahedron(1.0,5);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glColor4f(0.5f,0.5f,0.1f,0.67f);
	glDrawSphereIcosahedron(1.0,5);
	glDepthMask(GL_TRUE);
	glDisable(GL_BLEND);
	glEndList();
	}

/*************
Main function:
*************/

VRUI_APPLICATION_RUN(MDVisualizer)
