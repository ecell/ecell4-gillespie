#ifndef __GILLESPIE_SIMULATOR_HPP
#define __GILLESPIE_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <core/RandomNumberGenerator.hpp>
#include <core/Model.hpp>
#include <core/Simulator.hpp>

#include "GillespieWorld.hpp"
#include "GillespieSolver.hpp"

namespace ecell4 {

namespace gillespie {

class GillespieSimulator : public Simulator {
	/* XXX Memo */
	/*	Simulator class has t, num_steps, step, step(Real const &upto) funcs. */

public:
	GillespieSimulator(
		boost::shared_ptr<NetworkModel> model, boost::shared_ptr<GillespieWorld> world,
		RandomNumberGenerator &rng)
		: model_(model), world_(world), rng_(rng)
	{
		this->num_steps_ = 0;
	}
		
	Integer num_steps(void) const;
	void step(void) ;
	bool step(Real const & upto);
	void run(void);

	Real t(void) const;
	void set_t(Real const &t);

	RandomNumberGenerator &rng(void);
					
protected:
	boost::shared_ptr<NetworkModel> model_;
	boost::shared_ptr<GillespieWorld> world_;
	
	Integer num_steps_;
	RandomNumberGenerator &rng_;
};

}

}	// ecell4

#endif __GILLESPIE_SIMULATOR_HPP
