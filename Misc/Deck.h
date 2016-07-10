#ifndef _H_DECK_H_
#define _H_DECK_H_

#include "Card.h"
#include <list>
#include <algorithm>
#include <vector>

class Deck {
	/*	Class that will be used to fill the dealer's 'shoe' (multiple decks shuffled together to form a pile of cards that will be dealt)
	As many casinos and tournaments use either 6 or 8 decks in blackjack, this was generalised to form an n deck list based on player preference.

	The shoe is implemented as a list of type <Card>- this was chosen as the deck will constantly need to be resized as cards are dealt and
	removed from the pile, which would cause an overhead for other stl containers such as vectors. In addition, as only top card (ie the
	first card in the list) will need to be retrieved, the random access benefits of vectors will not be applicable here.

	As the 'shuffle' algorith requires random access, a vector of pairs<int card, int suit> will be created as a template for the shuffled
	card order.
	*/
private:
	std::list<Card> deck;
	int nDecks;
public:
	Deck(int _nDecks);
	Deck();
	Card drawCard();
	void newDeck();
	std::vector<std::pair<int, int>> shuffle();
	int getDeckSize();
};

#endif