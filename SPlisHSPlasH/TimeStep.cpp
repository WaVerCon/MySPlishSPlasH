#include "TimeStep.h"
#include "TimeManager.h"
#include "SPHKernels.h"
#include "PBF/SimulationDataPBF.h"
#include "SPlisHSPlasH/Utilities/Timing.h"
#include "SurfaceTension/SurfaceTension_Becker2007.h"
#include "SurfaceTension/SurfaceTension_Akinci2013.h"
#include "SurfaceTension/SurfaceTension_He2014.h"
#include "Viscosity/Viscosity_XSPH.h"
#include "Viscosity/Viscosity_Standard.h"


using namespace SPH;
using namespace std;

TimeStep::TimeStep(FluidModel *model)
{
	m_model = model;
	m_iterations = 0;
	m_cflMethod = 1;
	m_cflFactor = 0.5;
	m_cflMaxTimeStepSize = 0.005;
	m_maxIterations = 100;
	m_maxError = 0.01;
	m_maxIterationsV = 100;
	m_maxErrorV = 0.1;
	m_viscosity = NULL;
	setViscosityMethod(ViscosityMethods::XSPH); 
	m_surfaceTension = NULL;
	setSurfaceTensionMethod(SurfaceTensionMethods::None);
}

TimeStep::~TimeStep(void)
{
}

void TimeStep::clearAccelerations()
{
	const unsigned int count = m_model->numParticles();
	const Vector3r &grav = m_model->getGravitation();
	for (unsigned int i=0; i < count; i++)
	{
		// Clear accelerations of dynamic particles
		if (m_model->getMass(i) != 0.0)
		{
			Vector3r &a = m_model->getAcceleration(i);
			a = grav;
		}
	}
}

void TimeStep::computeDensities()
{
	const unsigned int numParticles = m_model->numParticles();
	const Real density0 = m_model->getDensity0();
	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			Real &density = m_model->getDensity(i);

			// Compute current density for particle i
			density = m_model->getMass(i) * m_model->W_zero();
			const Vector3r &xi = m_model->getPosition(0, i);

			for (unsigned int j = 0; j < m_model->numberOfNeighbors(i); j++)
			{
				const CompactNSearch::PointID &particleId = m_model->getNeighbor(i, j);
				const unsigned int &neighborIndex = particleId.point_id;
				const Vector3r &xj = m_model->getPosition(particleId.point_set_id, neighborIndex);

				if (particleId.point_set_id == 0)		// Test if fluid particle
				{
					density += m_model->getMass(neighborIndex) * m_model->W(xi - xj);
				}
				else 
				{
					// Boundary: Akinci2012
					density += m_model->getBoundaryPsi(particleId.point_set_id, neighborIndex) * m_model->W(xi - xj);
				}
			}
		}
	}
}

void TimeStep::updateTimeStepSize()
{
	if (m_cflMethod == 1)
		updateTimeStepSizeCFL(0.0001);
	else if(m_cflMethod == 2)
	{
		Real h = TimeManager::getCurrent()->getTimeStepSize();
		updateTimeStepSizeCFL(0.0001);
		if (m_iterations > 10)
			h *= 0.9;
		else if (m_iterations < 5)
			h *= 1.1;
		h = min(h, TimeManager::getCurrent()->getTimeStepSize());
		TimeManager::getCurrent()->setTimeStepSize(h);
	}
}

void TimeStep::updateTimeStepSizeCFL(const Real minTimeStepSize)
{
	const Real radius = m_model->getParticleRadius();
	Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Approximate max. position change due to current velocities
	Real maxVel = 0.1;
	const unsigned int numParticles = m_model->numParticles();
	const Real diameter = 2.0*radius;
	for (unsigned int i = 0; i < numParticles; i++)
	{
		const Vector3r &vel = m_model->getVelocity(0, i);
		const Vector3r &accel = m_model->getAcceleration(i);
		const Real velMag = (vel + accel*h).squaredNorm();
		if (velMag > maxVel)
			maxVel = velMag;
	}

	// boundary particles
	for (unsigned int i = 0; i < m_model->numberOfRigidBodyParticleObjects(); i++)
	{
		FluidModel::RigidBodyParticleObject *rbpo = m_model->getRigidBodyParticleObject(i);
		if (rbpo->m_rigidBody->isDynamic())
		{
			for (unsigned int j = 0; j < rbpo->numberOfParticles(); j++)
			{
				const Vector3r &vel = rbpo->m_v[j];
				const Real velMag = vel.squaredNorm();
				if (velMag > maxVel)
					maxVel = velMag;
			}
		}
	}

	// Approximate max. time step size 		
	h = m_cflFactor * .4 * (diameter / (sqrt(maxVel)));

	h = min(h, m_cflMaxTimeStepSize);
	h = max(h, minTimeStepSize);

	TimeManager::getCurrent()->setTimeStepSize(h);
}

void TimeStep::computeSurfaceTension()
{
	if (m_surfaceTension)
		m_surfaceTension->step();
}

void TimeStep::computeViscosity()
{
	if (m_viscosity)
		m_viscosity->step();
}

void TimeStep::performNeighborhoodSearch()
{
	START_TIMING("neighborhood_search");
	m_model->getNeighborhoodSearch()->find_neighbors();
	STOP_TIMING_AVG;
}

void TimeStep::reset()
{
	m_model->reset();
	if (m_surfaceTension)
		m_surfaceTension->reset();
	if (m_viscosity)
		m_viscosity->reset();
	m_iterations = 0;
}

void TimeStep::setSurfaceTensionMethod(SurfaceTensionMethods val)
{
	if ((val < SurfaceTensionMethods::None) || (val > SurfaceTensionMethods::He2014))
		val = SurfaceTensionMethods::None;

	if (val == m_surfaceTensionMethod)
		return;

	delete m_surfaceTension;
	m_surfaceTension = NULL;

	m_surfaceTensionMethod = val;
	if (m_surfaceTensionMethod == SurfaceTensionMethods::Becker2007)
		m_surfaceTension = new SurfaceTension_Becker2007(m_model);
	else if (m_surfaceTensionMethod == SurfaceTensionMethods::Akinci2013)
		m_surfaceTension = new SurfaceTension_Akinci2013(m_model);
	else if (m_surfaceTensionMethod == SurfaceTensionMethods::He2014)
		m_surfaceTension = new SurfaceTension_He2014(m_model);
}

void SPH::TimeStep::setViscosityMethod(ViscosityMethods val)
{
	if ((val < ViscosityMethods::None) || (val > ViscosityMethods::XSPH))
		val = ViscosityMethods::XSPH;

	if (val == m_viscosityMethod)
		return;

	delete m_viscosity;
	m_viscosity = NULL;

	m_viscosityMethod = val;

	if (m_viscosityMethod == ViscosityMethods::Standard)
		m_viscosity = new Viscosity_Standard(m_model);	
	else if (m_viscosityMethod == ViscosityMethods::XSPH)
		m_viscosity = new Viscosity_XSPH(m_model);
}

