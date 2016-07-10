#ifndef _H_CARD_H_
#define _H_CARD_H_

#include <string>
#include <vector>
#include <iostream>

class Card {
	// generic card class that will be used to populate the deck. Constructed with arguments (int card, int suit).
	// card and suit are zero based values (ie card = 0 represents an ace, card = 6 represents 7).
	// suit order is hearts, diamonds, clubs, spades
private:
	int card;
	int suit;
public:
	Card(int _card, int _suit);
	int getSuit();
	int getCard();
	void printCard();
	static const std::vector<int> CARDS_VAL_ARRAY;
	static const std::vector<std::string> CARDS_STR_ARRAY;
	static const std::vector<std::string> SUIT_STR_ARRAY;

};



#endif