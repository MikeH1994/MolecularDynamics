#include "GenericPlayer.h"

GenericPlayer::GenericPlayer(float money, std::string name) {
	this->money = money;
	this->score = 0;
	this->name = name;
	this->hasFolded = false;
	this->hasLost = false;
}

void GenericPlayer::clearHand() {
	hand.clear();
}

void GenericPlayer::addCard(Card card) {
	hand.push_back(card);
}

void GenericPlayer::addMoney(float winnings) {
	money += winnings;
	if (money < 0) {
		money = 0;
		hasLost = true;
	}
}

void GenericPlayer::printMove(Moves moves) {
	std::string str = getStrMoves(moves);
	std::cout << getName() << " " << str << std::endl;
}

void GenericPlayer::printHand() {
	std::cout << "Player " << getName() << " hand: ";
	for (unsigned int i = 0; i < hand.size(); i++) {
		std::cout << Card::CARDS_STR_ARRAY[hand[i].getCard()] << "; ";
	}
	std::cout << " score: " << getScore() << std::endl;
}

bool GenericPlayer::hasBlackjack() {
	if (hand.size() == 2 && getScore() == 21) {
		return true;
	}
	else return false;
}

std::string GenericPlayer::getName() {
	return this->name;
}

void GenericPlayer::setBet(float bet) {
	currentBet = bet;
	money -= bet;
}

void GenericPlayer::doubleDown(bool print) {
	money -= currentBet;
	currentBet *= 2;
	if (print) {
		std::cout << getName() << " doubles; current bet: " << currentBet << "; current money: " << getMoney() << std::endl;
	}
}

float GenericPlayer::getBet() {
	return currentBet;
}

std::vector<Card>& GenericPlayer::getHand() {
	return hand;
}

int GenericPlayer::getScore() {
	int score = 0;
	int nAcesHigh = 0; // number of aces the player has chosen to be 11

	for (unsigned int i = 0; i < hand.size(); i++) {
		Card card = hand[i];
		// add up the score using the default values (ie ace = 11; if it would cause player
		// to go bust, change this to one later)
		score += Card::CARDS_VAL_ARRAY[card.getCard()];
		// find number of aces
		if (card.getCard() == 0) {
			nAcesHigh += 1;
		}
	}
	if (score>21) {
		while (nAcesHigh > 0 && score>21) {
			score -= 10;
			nAcesHigh -= 1;
		}
	}

	return score;
}

std::string GenericPlayer::getStrMoves(Moves moves) {
	switch (moves) {
	case HIT:
		return "hits";
	case STICK:
		return "sticks";
	case DOUBLE:
		return "doubles";
	case QUIT:
		return "quits";
	case FOLD:
		return "folds";
	case BUST:
		return "has gone bust";
	}
	return "ERR";
}

float GenericPlayer::getMoney() {
	return money;
}