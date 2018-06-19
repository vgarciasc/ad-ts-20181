#include <iostream>
#include <queue>
#include "Event.h"
//#include <vld.h>

using namespace std;
const int TD = 500; // V.A. exp. lambda1 para tempo até a próxima chegada de dados

const int TVA = 16; // Tempo até o próximo pacote de voz durante o período ativo
const int TVS = 650; // V.A. exp. com média 650 para fim do período de silêncio do canal de voz
const double voiceContinue = 1.0 / 22; // Probabilidade de continuar o período ativo do canal de voz

const int AMOSTRAS = 1000;

int main(int argc, char *argv[]) {
	priority_queue<Event> arrivals;
	queue<Event> data, voice;
	// Calcula o tempo de chegada de cada evento
	arrivals.push(Event(TD, DATA));
	for (int i = 0; i < 30; ++i) {
		arrivals.push(Event(TVS, VOICE));
	}
	Event server(INTMAX_MAX, SERVER);
	for (int i = 0; i < AMOSTRAS; ++i) {
		Event arrival = arrivals.top();
		arrivals.pop();
		while (server < arrival) {
			// Freguês sai do servidor e puxa o próximo da fila
			if (!voice.empty()) {

			} else if (!data.empty()) {

			} else {

			}
		}

		// Adiciona a nova chegada na fila
		switch (arrival.type) {
			case DATA:
				data.push(arrival);
				arrivals.push(Event(arrival.time + TD, DATA));
				break;
			case VOICE:
				voice.push(arrival);
				if (voiceContinue)
					arrivals.push(Event(arrival.time + TVA, VOICE));
				else
					arrivals.push(Event(arrival.time + TVA + TVS, VOICE));
				break;
			default:
				break;
		}
	}
	return 0;
}