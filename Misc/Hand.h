#ifndef _H_HAND_H_
#define  _H_HAND_H_
#include "Card.h"


class Hand {
protected:
	std::vector<Card> hand;
public:
	void clearHand() {
		hand.clear();

	}

	void addCard(Card card) {
		hand.push_back(card);
	}

	bool hasBlackjack() {
		if (hand.size() == 2 && getScore() == 21) {
			return true;
		}
		else return false;
	}

};

#endif