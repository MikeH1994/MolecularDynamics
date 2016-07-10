#ifndef _H_HILO_H_
#define _H_HILO_H_

#include "Dealer.h"
#include "Card.h"
#include "Helper.h"

class HiLo {
	
public:
	HiLo(int nDecks);
	~HiLo();
protected:
	const enum Outcome {
		CORRECT, INCORRECT, DRAW
	};
	Deck* deck;
	int nCorrect = 0;
	int nIncorrect = 0;
	Outcome getOutcome(Card lastCard, Card nextCard, Helper::playerInput guessInput);
	void playHiLo();
	void updateCounter(Outcome outcome);

};

#endif