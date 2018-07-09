/*
 *
 */

#ifndef AD_STATS_H
#define AD_STATS_H

enum class PackageType {
	DATA, VOICE
};

struct ChannelProperties {
	int number;
	bool lastVoicePackage;
};

union Property {
	double wastedTime = 0;
	ChannelProperties channel;
};

class Packet {
public:
	int simulation = 0;
	double totalTime = 0;
	double serviceTime = 0;
	PackageType type;
	Property property;

	/**
	 * Construtor para pacotes de dados
	 * @param simulation
	 * @param serviceTime
	 */
	Packet(int simulation, double serviceTime);

	/**
	 * Construtor para pacotes de voz
	 * @param simulation
	 * @param channel
	 * @param serviceTime
	 */
	Packet(int simulation, int channel, double serviceTime);

	bool operator<(Packet other);
};


#endif
