/***********************************************************************
MDVisualizer - Vrui application class to show time-varying MD trace
data.
Copyright (c) 2005-2021 Oliver Kreylos
***********************************************************************/

#ifndef MDVISUALIZER_INCLUDED
#define MDVISUALIZER_INCLUDED

#define USE_IMPOSTORS 1

#include <utility>
#include <vector>
#include <Misc/StandardHashFunction.h>
#include <Misc/HashTable.h>
#include <Threads/MutexCond.h>
#include <Threads/Thread.h>
#include <Threads/Barrier.h>
#include <Threads/TripleBuffer.h>
#include <Geometry/ValuedPoint.h>
#include <Geometry/ArrayKdTree.h>
#include <GL/gl.h>
#include <GL/GLMaterial.h>
#include <GL/GLObject.h>
#if USE_IMPOSTORS
#include <GL/GLSphereRenderer.h>
#include <GL/GLCylinderRenderer.h>
#endif
#include <GL/GLGeometryVertex.h>
#include <GLMotif/ToggleButton.h>
#include <GLMotif/TextFieldSlider.h>
#include <Vrui/Application.h>
#include <Collaboration2/DataType.h>
#include <Collaboration2/Plugins/KoinoniaClient.h>

#include "MDGeometry.h"
#include "Atom.h"

/* Forward declarations: */
namespace Misc {
class CallbackData;
}
namespace Cluster {
class MulticastPipe;
}
namespace GLMotif {
class PopupMenu;
class PopupWindow;
}
namespace MD {
class Trace;
}

/* Namespace shortcuts: */
using Collab::DataType;
using Collab::Plugins::KoinoniaProtocol;
using Collab::Plugins::KoinoniaClient;

class MDVisualizer:public Vrui::Application,public GLObject
	{
	/* Embedded classes: */
	private:
	typedef Geometry::ValuedPoint<MD::Point,std::pair<unsigned int,MD::Scalar> > AtomPoint; // Value type is (atom index, half bond length)
	typedef Geometry::ArrayKdTree<AtomPoint> AtomTree;
	typedef std::pair<unsigned int,unsigned int> Bond; // Type for atom bonds
	typedef std::vector<Bond> BondList; // Type for lists of atom bonds
	
	class BondFinder; // Helper class to find bonds between atoms
	struct BondFinderThread; // Structure containing state for a worker thread to find bonds in loaded trace time steps
	
	struct TimeStep // Structure containing all data related to a single MD trace data time step
		{
		/* Elements: */
		public:
		size_t index; // Index of the time step
		size_t numAtoms; // Number of atoms in the time step
		MD::Atom* atoms; // Array of atoms (type ID and position)
		AtomTree atomTree; // Kd-tree storing all atom positions and indices for bond detection and selection
		BondList bonds; // Bonds between atoms (pairs of atom indices)
		
		/* Constructors and destructors: */
		TimeStep(void)
			:numAtoms(0),atoms(0)
			{
			};
		~TimeStep(void)
			{
			delete[] atoms;
			};
		};
	
	struct RenderingParameters // Structure for rendering parameters
		{
		/* Elements: */
		public:
		bool drawAtoms; // Flag whether to draw atoms
		float atomRadius; // Sphere radius to draw atoms
		bool drawBonds; // Flag whether to draw bonds
		float bondRadius; // Cylinder radius to draw bonds
		float selectedAtomRadius; // Sphere radius to draw selected atoms
		};
	
	struct DataItem:public GLObject::DataItem
		{
		/* Embedded classes: */
		public:
		typedef GLGeometry::Vertex<void,0,GLubyte,4,void,GLfloat,3> Vertex; // Type for rendered points
		
		/* Elements: */
		GLuint atomVertexBufferObjectId; // ID of vertex buffer object holding atom positions and colors
		GLuint bondIndexBufferObjectId; // ID of index buffer object holding bond atom pair indices
		unsigned int bufferVersion; // Version number of the data cached in the buffer objects
		GLuint influenceSphereDisplayListId; // ID of display list to render transparent spheres
		
		/* Constructors and destructors: */
		DataItem(void);
		virtual ~DataItem(void);
		};
	
	/* Elements: */
	private:
	MD::Trace* trace; // Visualized MD simulation trace
	const MD::Atom::TypeMap* atomTypeMap; // The trace's atom type map
	MD::Box boundingBox; // Bounding box of the entire simulation trace
	MD::Scalar maxHalfBondLength; // Maximum bond length defined for the entire simulation trace
	KoinoniaClient* koinonia; // Koinonia plug-in protocol client
	
	/* Background time step loader state: */
	Threads::MutexCond loadRequestCond; // A condition variable to request another time step from the loader thread
	bool loadRequest; // Flag if a load request has been issued
	size_t loadTimeStepIndex; // Index of the next time step to be loaded
	KoinoniaProtocol::ObjectID loadTimeStepIndexId; // Sharing ID of the index of the time stamp to load
	Cluster::MulticastPipe* timeStepPipe; // Pipe to stream trace time step data across a network in cluster rendering mode
	bool loaderThreadRunning; // Flag to keep the loader thread running
	Threads::Thread loaderThread; // Thread to load trace time steps in the background
	Threads::MutexCond bondFinderRequestCond; // Condition variable to wake up the bond finder threads
	bool bondFinderThreadsRunning; // Flag to keep the bond finder threads running
	unsigned int numBondFinderThreads; // Number of worker threads running in parallel to find bonds in a loaded time step
	BondFinderThread* bondFinderThreads; // Array of bond finder thread states
	Threads::Barrier bondFindersDone; // Barrier to wait for completion of all bond finder threads
	Threads::TripleBuffer<TimeStep> timeSteps; // Triple buffer of MD simulation trace time steps
	unsigned int timeStepVersion; // Version number of the currently locked trace buffer
	bool haveTimeStep; // Flag if a valid trace buffer is currently locked
	
	/* Interaction and visualization state: */
	Misc::HashTable<unsigned int,void> selectedAtomsSet; // Hash table containing indices of selected atoms
	std::vector<unsigned int> selectedAtoms; // List of indices of selected atoms
	GLMaterial atomMaterial; // Material to render selected atoms
	RenderingParameters renderingParameters; // The current rendering parameters
	KoinoniaProtocol::ObjectID renderingParametersId; // Sharing ID of the rendering parameters structure
	#if USE_IMPOSTORS
	GLSphereRenderer sphereRenderer; // Helper object to render atoms as spheres
	GLCylinderRenderer cylinderRenderer; // Helper object to render bonds as cylinders
	GLSphereRenderer selectedSphereRenderer; // Helper object to render selected atoms as spheres
	#endif
	
	double animationSpeed; // Speed of animation (in time steps per second)
	bool animate; // Flag whether to update trace files automatically
	double nextTimeStepTime; // Application time at which to display the next time step (if it has arrived)
	
	GLMotif::PopupMenu* mainMenu; // The program's main menu
	GLMotif::PopupWindow* animationDialog; // The animation control dialog
	GLMotif::TextFieldSlider* timeStepIndexSlider;
	GLMotif::ToggleButton* animateToggle;
	GLMotif::PopupWindow* renderDialog; // The rendering control dialog
	
	/* Private methods: */
	void* bondFinderThreadMethod(BondFinderThread* bondFinderThread); // Thread method to find bonds in loaded trace time steps in the background
	void* loaderThreadMethod(void); // Thread method to load trace time steps in the background
	void loadTimeStep(size_t newTimeStepIndex); // Requests loading the given time step
	static void loadTimeStepIndexUpdatedCallback(KoinoniaClient* client,KoinoniaProtocol::ObjectID id,void* object,void* userData);
	static void renderingParametersUpdatedCallback(KoinoniaClient* client,KoinoniaProtocol::ObjectID id,void* object,void* userData);
	void clearSelectionCallback(Misc::CallbackData* cbData);
	void showAnimationDialogCallback(Misc::CallbackData* cbData);
	void showRenderDialogCallback(Misc::CallbackData* cbData);
	GLMotif::PopupMenu* createMainMenu(void);
	void timeStepIndexCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData);
	void animationSpeedCallback(GLMotif::TextFieldSlider::ValueChangedCallbackData* cbData);
	void animateCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	GLMotif::PopupWindow* createAnimationDialog(void);
	void renderingParametersChangedCallback(Misc::CallbackData* cbData);
	GLMotif::PopupWindow* createRenderDialog(void);
	
	/* Constructors and destructors: */
	public:
	MDVisualizer(int& argc,char**& argv);
	virtual ~MDVisualizer(void);
	
	/* Methods from class Vrui::Application: */
	virtual void frame(void);
	virtual void display(GLContextData& contextData) const;
	virtual void resetNavigation(void);
	
	/* Methods from class GLObject: */
	virtual void initContext(GLContextData& contextData) const;
	};

#endif
