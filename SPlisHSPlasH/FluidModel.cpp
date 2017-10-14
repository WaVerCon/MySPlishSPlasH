#include "FluidModel.h"
#include "SPHKernels.h"
#include <iostream>

using namespace SPH;

FluidModel::FluidModel() 
{	
	m_density0 = 1000.0;
	m_particleRadius = 0.025;
	m_viscosity = 0.02;
	m_neighborhoodSearch = NULL;
	m_gravitation = Vector3r(0.0, -9.81, 0.0);
	m_stiffness = 50000.0;
	m_exponent = 7.0;
	m_surfaceTension = 0.05;
	m_enableDivergenceSolver = true;
	m_velocityUpdateMethod = 0;

	ParticleObject *fluidParticles = new ParticleObject();
	m_particleObjects.push_back(fluidParticles);

	setKernel(0);
	setGradKernel(0);
}

FluidModel::~FluidModel(void)
{
	cleanupModel();
}

void FluidModel::cleanupModel()
{
	releaseFluidParticles();
	for (unsigned int i = 0; i < m_particleObjects.size(); i++)
	{
		if (i > 0)
		{
			RigidBodyParticleObject *rbpo = ((RigidBodyParticleObject*)m_particleObjects[i]);
			rbpo->m_boundaryPsi.clear();
			rbpo->m_f.clear();
			delete rbpo->m_rigidBody;
			delete rbpo;
		}
		else
			delete m_particleObjects[i];
	}
	m_particleObjects.clear();

	m_a.clear();
	m_masses.clear();
	m_density.clear();
	delete m_neighborhoodSearch;
}

void FluidModel::reset()
{
	const unsigned int nPoints = numParticles();

	// reset velocities and accelerations
	for (unsigned int i = 1; i < m_particleObjects.size(); i++)
	{
		for (int j = 0; j < (int)m_particleObjects[i]->m_x.size(); j++)
		{
			RigidBodyParticleObject *rbpo = ((RigidBodyParticleObject*)m_particleObjects[i]);
			rbpo->m_f[j].setZero();
			rbpo->m_v[j].setZero();
		}
	}
	
	// Fluid
	for(unsigned int i=0; i < nPoints; i++)
	{
		const Vector3r& x0 = getPosition0(0, i);
		getPosition(0, i) = x0;
		getVelocity(0, i).setZero();
		getAcceleration(i).setZero();
		m_density[i] = 0.0;
	}

	if (m_neighborhoodSearch != NULL)
	{
		performNeighborhoodSearchSort();
	}
	updateBoundaryPsi();
}

void FluidModel::initMasses()
{
	const int nParticles = (int) numParticles();
	const Real diam = 2.0*m_particleRadius;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < nParticles; i++)
		{
			setMass(i, 0.8 * diam*diam*diam * m_density0);		// each particle represents a cube with a side length of r		
																// mass is slightly reduced to prevent pressure at the beginning of the simulation
		}
	}
}

void FluidModel::resizeFluidParticles(const unsigned int newSize)
{
	m_particleObjects[0]->m_x0.resize(newSize);
	m_particleObjects[0]->m_x.resize(newSize);
	m_particleObjects[0]->m_v.resize(newSize);
	m_a.resize(newSize);
	m_masses.resize(newSize);
	m_density.resize(newSize);
}

void FluidModel::releaseFluidParticles()
{
	m_particleObjects[0]->m_x0.clear();
	m_particleObjects[0]->m_x.clear();
	m_particleObjects[0]->m_v.clear();
	m_a.clear();
	m_masses.clear();
	m_density.clear();
}

void FluidModel::initModel(const unsigned int nFluidParticles, Vector3r* fluidParticles)
{
	releaseFluidParticles();
	resizeFluidParticles(nFluidParticles);

	// init kernel
	setParticleRadius(m_particleRadius);

	// copy fluid positions
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)nFluidParticles; i++)
		{
			getPosition0(0, i) = fluidParticles[i];
		}
	}

	// initialize masses
	initMasses();

	// Initialize neighborhood search
	if (m_neighborhoodSearch == NULL)
		m_neighborhoodSearch = new CompactNSearch::NeighborhoodSearch(m_supportRadius, false);
	m_neighborhoodSearch->set_radius(m_supportRadius);

	// Fluids 
	m_neighborhoodSearch->add_point_set(&getPosition(0, 0)[0], nFluidParticles, true, true);

	// Boundary
	for (unsigned int i = 0; i < numberOfRigidBodyParticleObjects(); i++)
	{
		RigidBodyParticleObject *rb = getRigidBodyParticleObject(i);
		m_neighborhoodSearch->add_point_set(&rb->m_x[0][0], rb->m_x.size(), rb->m_rigidBody->isDynamic(), false);
	}

	reset();
}
/*
在初始化流体模型，改变当前流体模型和构造TimeStep时将调用该函数。
先激活静态边界并计算边界粒子系数boundaryPsi，再分别激活每个动态边界并计算边界粒子系数。最后只激活流体粒子，非激活边界粒子。
激活指是否参与领域搜索。
参考NeighborhoodSearch
*/
void FluidModel::updateBoundaryPsi()
{
	//////////////////////////////////////////////////////////////////////////
	// Compute value psi for boundary particles (boundary handling)
	// (see Akinci et al. "Versatile rigid - fluid coupling for incompressible SPH", Siggraph 2012
	//////////////////////////////////////////////////////////////////////////

	// Search boundary neighborhood

	// Activate only static boundaries
	std::cout << "Initialize boundary psi\n";
	m_neighborhoodSearch->point_set(0).enable_neighborsearch(false);
	for (unsigned int i = 0; i < numberOfRigidBodyParticleObjects(); i++)
	{
		if (!getRigidBodyParticleObject(i)->m_rigidBody->isDynamic())
			m_neighborhoodSearch->point_set(i + 1).enable_neighborsearch(true);
	}

	m_neighborhoodSearch->find_neighbors();

	// Boundary objects
	for (unsigned int body = 0; body < numberOfRigidBodyParticleObjects(); body++)
	{
		if (!getRigidBodyParticleObject(body)->m_rigidBody->isDynamic())
			computeBoundaryPsi(body);
	}

	////////////////////////////////////////////////////////////////////////// 
	// Compute boundary psi for all dynamic bodies
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int body = 0; body < numberOfRigidBodyParticleObjects(); body++)
	{
		// Deactivate all
		for (int j = 1; j < m_neighborhoodSearch->point_sets().size(); j++)
			m_neighborhoodSearch->point_set(j).enable_neighborsearch(false);

		// Only activate next dynamic body
		if (getRigidBodyParticleObject(body)->m_rigidBody->isDynamic())
		{
			m_neighborhoodSearch->point_set(body + 1).enable_neighborsearch(true);
			m_neighborhoodSearch->find_neighbors();
			computeBoundaryPsi(body);
		}
	}

	// Activate only fluids
	m_neighborhoodSearch->point_set(0).enable_neighborsearch(true);
	for (int i = 1; i < m_neighborhoodSearch->point_sets().size(); i++)
		m_neighborhoodSearch->point_set(i).enable_neighborsearch(false);

}
/*
对边界粒子搜索领域粒子，根据领域粒子密度计算boundaryPsi。
关键方法：NeighborhoodSearch->point_set.neighbor(i,j)，在指定的点集（按照先前ParticleObject进行添加，每个PartcileObject是一个点集）内搜索邻居。
计算边界粒子对流体粒子的贡献。详见Akinci[2012]
*/
void FluidModel::computeBoundaryPsi(const unsigned int body)
{
	const Real density0 = getDensity0();
	
	RigidBodyParticleObject *rb = getRigidBodyParticleObject(body);
	const unsigned int numBoundaryParticles = rb->numberOfParticles();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numBoundaryParticles; i++)
		{
			Real delta = m_W_zero;
			for (unsigned int j = 0; j < m_neighborhoodSearch->point_set(body + 1).n_neighbors(i); j++)
			{
				const CompactNSearch::PointID &pid = m_neighborhoodSearch->point_set(body + 1).neighbor(i, j);
				if (pid.point_set_id != 0)
					delta += W(getPosition(body + 1, i) - getPosition(pid.point_set_id, pid.point_id));
			}
			const Real volume = 1.0 / delta;
			rb->m_boundaryPsi[i] = density0 * volume; 
		}
	}
}

void FluidModel::addRigidBodyObject(RigidBodyObject *rbo, const unsigned int numBoundaryParticles, Vector3r *boundaryParticles)
{
	//RigidBodyParticleObject对象内部有数个边界粒子，并包含位置、静止位置、速度、力。
	RigidBodyParticleObject *rb = new RigidBodyParticleObject();
	m_particleObjects.push_back(rb);

	rb->m_x0.resize(numBoundaryParticles);
	rb->m_x.resize(numBoundaryParticles);
	rb->m_v.resize(numBoundaryParticles);
	rb->m_f.resize(numBoundaryParticles);
	rb->m_boundaryPsi.resize(numBoundaryParticles);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numBoundaryParticles; i++)
		{
			rb->m_x0[i] = boundaryParticles[i];
			rb->m_x[i] = boundaryParticles[i];
			rb->m_v[i].setZero();
			rb->m_f[i].setZero();
		}
	}
	rb->m_rigidBody = rbo;//RigidBodyParticleObject对象有一个指向存储刚体对象属性的RigidBodyObject（子类有PBDRigidBody）指针。
}
//调用neighborhoodSearch索引粒子
void FluidModel::performNeighborhoodSearchSort()
{
	const unsigned int numPart = numParticles();
	if (numPart == 0)
		return;

	m_neighborhoodSearch->z_sort();

	auto const& d = m_neighborhoodSearch->point_set(0);
	d.sort_field(&m_particleObjects[0]->m_x0[0]);
	d.sort_field(&m_particleObjects[0]->m_x[0]);
	d.sort_field(&m_particleObjects[0]->m_v[0]);
	d.sort_field(&m_a[0]);
	d.sort_field(&m_masses[0]);
	d.sort_field(&m_density[0]);


	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int i = 1; i < m_neighborhoodSearch->point_sets().size(); i++)
	{
		RigidBodyParticleObject *rb = getRigidBodyParticleObject(i - 1);// m_particleOjbects[0]为流体粒子，m_particleOjbects[i](i>0)是刚体粒子
		if (rb->m_rigidBody->isDynamic())			// sort only dynamic boundaries
		{
			auto const& d = m_neighborhoodSearch->point_set(i);
			d.sort_field(&rb->m_x0[0]);
			d.sort_field(&rb->m_x[0]);
			d.sort_field(&rb->m_v[0]);
			d.sort_field(&rb->m_f[0]);
			d.sort_field(&rb->m_boundaryPsi[0]);
		}
	}
}

void SPH::FluidModel::setParticleRadius(Real val)
{
	m_particleRadius = val; 
	m_supportRadius = 4.0*m_particleRadius;

	// init kernel
	Poly6Kernel::setRadius(m_supportRadius);
	SpikyKernel::setRadius(m_supportRadius);
	CubicKernel::setRadius(m_supportRadius);
	PrecomputedCubicKernel::setRadius(m_supportRadius);
	CohesionKernel::setRadius(m_supportRadius);
	AdhesionKernel::setRadius(m_supportRadius);
}


void SPH::FluidModel::setGradKernel(unsigned int val)
{
	m_gradKernelMethod = val;
	if (m_gradKernelMethod == 0)
		m_gradKernelFct = CubicKernel::gradW;
	else if (m_gradKernelMethod == 1)
		m_gradKernelFct = Poly6Kernel::gradW;
	else if (m_gradKernelMethod == 2)
		m_gradKernelFct = SpikyKernel::gradW;
	else if (m_gradKernelMethod == 3)
		m_gradKernelFct = FluidModel::PrecomputedCubicKernel::gradW;
}

void SPH::FluidModel::setKernel(unsigned int val)
{
	m_kernelMethod = val;
	if (m_kernelMethod == 0)
	{
		m_W_zero = CubicKernel::W_zero();
		m_kernelFct = CubicKernel::W;
	}
	else if (m_kernelMethod == 1)
	{
		m_W_zero = Poly6Kernel::W_zero();
		m_kernelFct = Poly6Kernel::W;
	}
	else if (m_kernelMethod == 2)
	{
		m_W_zero = SpikyKernel::W_zero();
		m_kernelFct = SpikyKernel::W;
	}
	else if (m_kernelMethod == 3)
	{
		m_W_zero = FluidModel::PrecomputedCubicKernel::W_zero();
		m_kernelFct = FluidModel::PrecomputedCubicKernel::W;
	}
}