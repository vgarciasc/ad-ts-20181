#include "Event.h"

Event::Event(double time, EventType type, Stats* stats) {
	this->time = time;
	this->type = type;
	this->stats = stats;
}
Event::Event(double time, EventType type) {
	this->time = time;
	this->type = type;
	this->stats = new Stats();
}

//Event::Event(double time, double serviceTime, EventType type) {
//	this->time = time;
//	this->type = type;
//	this->serviceTime = serviceTime;
//}

bool
operator<(const Event lhs, const Event rhs) {
	return lhs.time < rhs.time;
}

//Event::Event(double time, double serviceTime, EventType type, int channel) {
//	this->time = time;
//	this->type = type;
//	this->serviceTime = serviceTime;
//	this->channel = channel;
//}
