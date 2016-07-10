#ifndef _H_GENERICPLAYER_H_
#define _H_GENERICPLAYER_H_
#include "Card.h"

#include <vector>
#include <string>

class GenericPlayer {
protected:
	float money = 0;
	float currentBet = 0;
	int score = 0;
	std::vector<Card> hand;
	bool hasFolded;
	bool hasLost;

	std::string name;

public:
	const enum Moves {
		FOLD, HIT, DOUBLE, STICK, SPLIT, QUIT, BUST, ERR
	};
	const enum PlayerType {
		NPC_SIDNEY, NPC_ALEX
	};

	GenericPlayer(float money, std::string name);
	GenericPlayer() {}
	void clearHand();
	void addCard(Card card);
	void addMoney(float winnings);
	void printMove(Moves moves);
	void printHand();
	bool hasBlackjack();
	std::string getName();
	void setBet(float bet);
	void doubleDown(bool print);
	float getBet();
	virtual Moves getMove(bool print) = 0;
	std::vector<Card>& getHand();
	int getScore();
	std::string getStrMoves(Moves moves);
	float getMoney();

};


#endif
