#include "Dealer.h"
#include "GenericPlayer.h"

Dealer::Dealer() {
	this->deck = Deck(1);
}

Dealer::Dealer(int nDecks) :GenericPlayer(1000, "Dealer") {
	this->deck = Deck(nDecks);
	this->nDecks = nDecks;
}

Card Dealer::drawNext() {
	return deck.drawCard();
}

GenericPlayer::Moves Dealer::getMove(bool print) {
	Moves move;
	//rules of dealer's moves; hit if score is below 7, stick if 
	if (getScore() < 17) {
		move = HIT;
	}
	else {
		move = STICK;
	}
	if (print) {
		printMove(move);
	}
	return move;
}

void Dealer::dealSelf(bool print){
	//Dealer deals self two cards- one of which is face down (hole card)
	Card card = drawNext();
	if (print){
		std::cout << "Dealer has dealt hole card" << std::endl;
	}
	addCard(card);
	card = drawNext();
	if (print) {
		std::cout << "Dealer has dealt ";
		card.printCard();
		std::cout << std::endl;
		std::cout << "------------------" << std::endl;
	}
	addCard(card);
}

template <class T>
void Dealer::dealCards(int nCards, T &player, bool print) {
	if (print) {
		std::cout << "Player " << player.getName() << " was dealt ";
	}
	for (int i = 0; i<nCards; i++) {
		Card card = drawNext();
		
		if (print) {
			card.printCard();
			std::cout << "; ";
		}
		player.addCard(card);
	}
	if (print) {
		std::cout << std::endl;
		std::cout << player.getName() << " score: " << player.getScore() << std::endl;
		if (player.getScore()>21) {
			std::cout << "======Bust======" << std::endl;
		}
		std::cout << "------------------" << std::endl;
	}
}

int Dealer::getDeckSize(){
	return deck.getDeckSize();
}

void Dealer::newDeck(){
	deck = Deck(nDecks);
}
