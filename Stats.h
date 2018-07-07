/*
 *
 */

#ifndef AD_STATS_H
#define AD_STATS_H


class Stats {
public:
	double enterQueueTime = 0;
	double waitTime = 0;
	double serviceTime = 0;
	double wastedServiceTime = 0;
	bool lastVoicePackage = false;
	int channel = 0;
	Stats();
	Stats(int channel, double serviceTime);
};


#endif
