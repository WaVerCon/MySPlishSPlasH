#ifndef __SimulationDataPCISPH_h__
#define __SimulationDataPCISPH_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Predictive-corrective Incompressible SPH introduced
	* by Solenthaler and Pajarola \cite Solenthaler:2009.
	*/
	class SimulationDataPCISPH
	{
		public:
			SimulationDataPCISPH();
			virtual ~SimulationDataPCISPH();

		protected:	
			FluidModel *m_model;
			Real m_pcisph_factor;
			//Matrix3r m_stress_factor;//stress scaling factor

			std::vector<Vector3r> m_lastX;
			std::vector<Vector3r> m_lastV;
			std::vector<Real> m_densityAdv;
			std::vector<Real> m_pressure;
			std::vector<Vector3r> m_pressureAccel;
			//std::vector<Matrix3r> m_stress;//record stress tensor

		public:
			/** Initialize the arrays containing the particle data.
			*/
			virtual void init(FluidModel *model);

			/** Release the arrays containing the particle data.
			*/
			virtual void cleanup();

			/** Reset the particle data.
			*/
			virtual void reset();

			/** Important: First call m_model->performNeighborhoodSearchSort() 
			 * to call the z_sort of the neighborhood search.
			 */
			void performNeighborhoodSearchSort();

			Real getPCISPH_ScalingFactor() { return m_pcisph_factor; }
			//Matrix3r getStressScalingFactor() { return m_stress_factor; }

			FORCE_INLINE Vector3r &getLastPosition(const unsigned int i)
			{
				return m_lastX[i];
			}

			FORCE_INLINE const Vector3r &getLastPosition(const unsigned int i) const
			{
				return m_lastX[i];
			}

			FORCE_INLINE void setLastPosition(const unsigned int i, const Vector3r &pos)
			{
				m_lastX[i] = pos;
			}

			FORCE_INLINE Vector3r &getLastVelocity(const unsigned int i)
			{
				return m_lastV[i];
			}

			FORCE_INLINE const Vector3r &getLastVelocity(const unsigned int i) const
			{
				return m_lastV[i];
			}

			FORCE_INLINE void setLastVelocity(const unsigned int i, const Vector3r &vel)
			{
				m_lastV[i] = vel;
			}

			FORCE_INLINE const Real getDensityAdv(const unsigned int i) const
			{
				return m_densityAdv[i];
			}

			FORCE_INLINE Real& getDensityAdv(const unsigned int i)
			{
				return m_densityAdv[i];
			}

			FORCE_INLINE void setDensityAdv(const unsigned int i, const Real d)
			{
				m_densityAdv[i] = d;
			}

			FORCE_INLINE const Real getPressure(const unsigned int i) const
			{
				return m_pressure[i];
			}

			FORCE_INLINE Real& getPressure(const unsigned int i)
			{
				return m_pressure[i];
			}

			FORCE_INLINE void setPressure(const unsigned int i, const Real p)
			{
				m_pressure[i] = p;
			}
			//stress getter and setter
			/*FORCE_INLINE const Matrix3r getStress(const unsigned int i) const
			{
				return m_stress[i];
			}

			FORCE_INLINE Matrix3r& getStress(const unsigned int i)
			{
				return m_stress[i];
			}

			FORCE_INLINE void setStress(const unsigned int i, const Matrix3r p)
			{
				m_stress[i] = p;
			}
			*/
			FORCE_INLINE Vector3r &getPressureAccel(const unsigned int i)
			{
				return m_pressureAccel[i];
			}

			FORCE_INLINE const Vector3r &getPressureAccel(const unsigned int i) const
			{
				return m_pressureAccel[i];
			}

			FORCE_INLINE void setPressureAccel(const unsigned int i, const Vector3r &val)
			{
				m_pressureAccel[i] = val;
			}

	};
}

#endif