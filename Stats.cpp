/*
 *
 */

#include "Stats.h"

Stats::Stats() {

}

Stats::Stats(int channel, double serviceTime) {
	this->channel = channel;
	this->serviceTime = serviceTime;
}
