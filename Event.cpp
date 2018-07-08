#include "Event.h"

Event::Event(double time, EventType type, Packet* stats) {
	this->time = time;
	this->type = type;
	this->packet = stats;
}

bool
operator<(const Event lhs, const Event rhs) {
	return rhs.time < lhs.time;
}
