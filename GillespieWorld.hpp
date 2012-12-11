#ifndef INCLUDE_GUARD_GILLESPIE_WORLD
#	define INCLUDE_GUARD_GILLESPIE_WORLD

#include <map>
#include <boost/scoped_ptr.hpp>
#include <string>

#include <core/CompartmentSpace.hpp>
#include <core/Species.hpp>


//============================================================
//	World	
//============================================================
class World {
private:
public:	// Data members
	double current_t;
	std::map<std::string,int> current_state;	// [id] -> number of substrate	// id => 分子の個数　

public:
	World();
	double get_current_time(void);
	void set_current_time(double new_t);
	int get_current_state(std::string &sp);
	void set_current_state(std::string &sp, int number);
	void add_specie(std::string &sp, int number);

	std::string to_string(void);
};

std::ostream &operator<<(std::ostream &s, World &w);

//============================================================
//	New Implementation
//============================================================
namespace ecell4 {
	
namespace gillespie {

class GillespieWorld {
public:
	GillespieWorld(Real const &volume)
		: cs_(new CompartmentSpaceVectorImpl(volume))
	{
		;
	}
	// about time
	void set_t(Real const &t);
	Real t(void);

	// about molecules states
	// 		immutable functions.
	Integer num_species(void);
	bool has_species(Species const &sp);
	Integer num_molecules(Species const& sp);
	
	//		mutable functions.
	void add_species(Species const &sp);
	void remove_species(Species const &sp);
	void add_molecules(Species const &sp, Integer const &num);
	// I think it is better that the name of this function is 'decrease_molecules()'.
	void remove_molecules(Species const &sp, Integer const &num);
	
private:
	boost::scoped_ptr<CompartmentSpace> cs_;
};
		

}	// gillespie

}	// ecell4


#endif	//INCLUDE_GUARD_GILLESPIE_WORLD
