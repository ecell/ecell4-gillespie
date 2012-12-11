#include <vector>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
#include "./GillespieWorld.hpp"

//============================================================
//	World	
//		-- General
//============================================================
World::World(void)
{
	this->current_t = 0.0;
}
double World::get_current_time(void)
{	return this->current_t;	}

void World::set_current_time(double new_t)
{	this->current_t = new_t;	}

int World::get_current_state(std::string &sp)
{	
	return this->current_state[sp];	
}

void World::set_current_state(std::string &sp, int number)
{	
	this->current_state[sp] = number;	
}

void World::add_specie(std::string &sp, int number = 0)
{
	this->current_state.insert(std::map<string,int>::value_type(sp, number));
}

//============================================================
//	World
//		-- pretty printer
//============================================================
string World::to_string(void) {
	ostringstream os;
	os << "time: " << this->current_t;
	for(std::map<string,int>::iterator it = this->current_state.begin(); it != this->current_state.end(); it++) {
		os << ", " << it->first << ": " << it->second;
	}
	return os.str();
}

ostream &operator<<(ostream &s, World &w) {
	return s << w.to_string();
}

//============================================================
//	New Implementation
//============================================================
namespace ecell4 {
	
namespace gillespie {
//============================================================
//			About Time Operation
//============================================================
void GillespieWorld::set_t(Real const &t) {
	this->cs_->set_t(t);
}
Real GillespieWorld::t(void) {
	return this->cs_->t();
}

//============================================================
//			About Molecules Operation
//============================================================
Integer GillespieWorld::num_species(void) {
	return this->cs_->num_species();
}

bool GillespieWorld::has_species(Species const &sp) {
	return this->cs_->has_species(sp);
}

Integer GillespieWorld::num_molecules(Species const &sp) {
	return this->cs_->num_molecules(sp);
}

void GillespieWorld::add_species(Species const &sp) {
	this->cs_->add_species(sp);
	return;
}

void GillespieWorld::remove_species(Species const &sp) {
	this->cs_->remove_species(sp);
	return;
}

void GillespieWorld::add_molecules(Species const &sp, Integer const &num) {
	this->cs_->add_molecules(sp, num);
	return;
}

void GillespieWorld::remove_molecules(Species const &sp, Integer const &num) {
	this->cs_->remove_molecules(sp, num);
	return;
}

}	//gillespie

}	//ecell4
