/*
 *
 */

#ifndef AD_EVENTO_H
#define AD_EVENTO_H
enum class EventType {
	EMPTY, SERVER, DATA, VOICE
};

class Event {
public:
	double time;
	double serviceTime;
	EventType type;
	Event(double time, EventType type);
	Event(double time, double serviceTime, EventType type);

	friend bool operator<(Event lhs, Event rhs);
};

#endif AD_EVENTO_H
