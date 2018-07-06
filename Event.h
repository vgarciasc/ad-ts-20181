/*
 *
 */

#ifndef AD_EVENTO_H
#define AD_EVENTO_H

#include "Stats.h"

enum class EventType {
	EMPTY, SERVER, DATA, VOICE
};

class Event {
public:
	double time = 0;
	EventType type;
	Stats* stats;
	Event(double time, EventType type);
	Event(double time, EventType type, Stats* stats);
//	Event(double time, double serviceTime, EventType type);
//	Event(double time, double serviceTime, EventType type, int channel);

	friend bool operator<(Event lhs, Event rhs);
};

#endif