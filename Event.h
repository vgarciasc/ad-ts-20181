/*
 *
 */

#ifndef AD_EVENTO_H
#define AD_EVENTO_H
enum EventType {
	SERVER, DATA, VOICE
};

class Event {
public:
	long long int time;
	EventType type;

	Event(long long int time, EventType type);

	friend bool operator<(Event lhs, Event rhs);
};

#endif AD_EVENTO_H
