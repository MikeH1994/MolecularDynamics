#include "Deck.h"
#include "Card.h"

Deck::Deck(int _nDecks):nDecks(_nDecks) {
	newDeck();
}

Deck::Deck() {
	nDecks = 1;
	newDeck();
}

void Deck::newDeck() {
	//using the shuffled order, generate the list 
	deck.clear();
	std::vector<std::pair<int, int>> shuffledOrder = shuffle();
	for (unsigned int i = 0; i < shuffledOrder.size(); i++) {
		std::pair<int, int> params = shuffledOrder[i];
		deck.push_back(Card(params.first, params.second));
	}
}

Card Deck::drawCard() {
	Card card = *deck.begin();
	deck.pop_front();
	return card;
}

std::vector<std::pair<int, int>> Deck::shuffle() {
	int nShuffles = 6;
	//creates a 'template' of how the shuffled card deck will look, this will be passed on in order to create a list of cards with this order.
	std::vector<std::pair<int, int>> cardOrder(52 * nDecks);
	int index = 0;//running count of position in vector
	for (int n = 0; n < nDecks; n++) {		// for every deck
		for (int j = 0; j < 4; j++) {		// for every suit
			for (int i = 0; i < 13; i++) {  // for every card
				cardOrder[index] = std::pair<int, int>(i, j);
				index++;
			}
		}
	}
	for (int i = 0; i < nShuffles; i++) {
		std::random_shuffle(cardOrder.begin(), cardOrder.end());
	}
	return cardOrder;
}

int Deck::getDeckSize() {
	return (int) deck.size();
}