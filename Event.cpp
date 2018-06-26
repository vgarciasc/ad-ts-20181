//
// Created by joao.albuquerque on 19/06/2018.
//

#include "Event.h"

Event::Event(
		double time,
		EventType type) {
	this->time = time;
	this->type = type;
}

Event::Event(
		double time,
		double serviceTime,
		EventType type) {
	this->time = time;
	this->type = type;
	this->serviceTime = serviceTime;
}

bool
operator<(
		const Event lhs,
		const Event rhs) {
	return lhs.time <
		   rhs.time;
}