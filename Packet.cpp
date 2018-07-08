/*
 *
 */

#include "Packet.h"

Packet::Packet(int simulation, double serviceTime) {
	this->simulation = simulation;
	this->serviceTime = serviceTime;
	this->type = PackageType::DATA;
}

Packet::Packet(int simulation, int channel, double serviceTime) {
	this->simulation = simulation;
	this->property.channel.number = channel;
	this->serviceTime = serviceTime;
	this->type = PackageType::VOICE;
}

bool Packet::operator<(Packet other) {
	if (this->type == other.type){
		return this->totalTime < other.totalTime;
	} else {
		return this->type == PackageType::VOICE;
	}
}