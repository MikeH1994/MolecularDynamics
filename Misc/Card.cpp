#include "Card.h"
#include <string>

Card::Card(int _card, int _suit) :card(_card), suit(_suit) {
	//printCard();
}

int Card::getSuit() {
	return this->suit;
}
int Card::getCard() {
	return this->card;
}
void Card::printCard() {
	//prints out "Ace of Spades", for example, based on the values of 'card' and 'suit'
	std::cout << CARDS_STR_ARRAY[card]<<" of "<<SUIT_STR_ARRAY[suit];
}

