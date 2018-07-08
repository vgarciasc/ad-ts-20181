/*
 *
 */

#ifndef AD_EVENTO_H
#define AD_EVENTO_H

#include "Packet.h"

enum class EventType {
	SERVER, DATA, VOICE
};

class Event {
public:
	double time = 0;
	EventType type;
	Packet* packet;
	Event(double time, EventType type, Packet* stats);

	friend bool operator<(Event lhs, Event rhs);
};

#endif