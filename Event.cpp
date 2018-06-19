//
// Created by joao.albuquerque on 19/06/2018.
//

#include "Event.h"

Event::Event(long long int time, EventType type) {
	this->time = time;
	this->type = type;
}

bool operator<(const Event lhs, const Event rhs) {
		return lhs.time < rhs.time;
}

