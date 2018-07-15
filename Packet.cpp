/*
 *
 */

#include "Packet.h"

Packet::Packet(int simulation, double serviceTime) {
	this->simulation = simulation;
	this->serviceTime = serviceTime;
	this->type = PacketType::DATA;
}

Packet::Packet(int simulation, int channel, double serviceTime) {
	this->simulation = simulation;
	this->property.channel.number = channel;
	this->serviceTime = serviceTime;
	this->type = PacketType::VOICE;
}

bool Packet::operator<(Packet other) {
	if (this->type == other.type){
		return this->totalTime < other.totalTime;
	} else {
		return this->type == PacketType::VOICE;
	}
}