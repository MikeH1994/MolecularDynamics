#include "HiLo.h"

HiLo::HiLo(int nDecks) {
	this->deck = new Deck(nDecks);
	playHiLo();
}

HiLo::~HiLo() {
	delete deck;
}

HiLo::Outcome HiLo::getOutcome(Card lastCard, Card nextCard, Helper::playerInput guessInput) {
	// returns 'correct', 'incorrect' or 'draw' based on cards and guess
	if (lastCard.getCard() == nextCard.getCard()) {
		return DRAW;
	}

	bool guess = false;
	if (guessInput == Helper::HI) {
		guess = true;
	}
	else {
		guess = false;
	}

	if (guess == (nextCard.getCard() > lastCard.getCard())) {
		return CORRECT;
	}
	else {
		return INCORRECT;
	}
}

void HiLo::playHiLo() {
	std::cout << "\n===================\nPlaying HiLo\n===================\n" << std::endl;
	std::vector<Helper::playerInput> allowedInput = { Helper::HI, Helper::LO, Helper::QUIT };
	Outcome outcome;
	Helper::playerInput guess = Helper::HI;

	Card lastCard = deck->drawCard();

	std::cout << "First card is ";
	lastCard.printCard();

	while (guess != Helper::QUIT) {

		Card nextCard = deck->drawCard();
		guess = Helper::getInput(allowedInput, "'HI' or 'LO'?");


		if (guess == Helper::QUIT) {
			break;
		}
		std::cout << "Next card is ";
		nextCard.printCard();
		outcome = getOutcome(lastCard, nextCard, guess);
		updateCounter(outcome);
		std::cout << "Correct: " << nCorrect << "; Incorrect: " << nIncorrect << std::endl;
		lastCard = nextCard;

		if (deck->getDeckSize() < 10) {
			deck = new Deck(2);
		}

	}
	std::cout << "End of Hi Lo game" << std::endl;
	this->~HiLo();
}

void HiLo::updateCounter(Outcome outcome) {
	if (outcome == CORRECT) {
		nCorrect++;
		std::cout << "Correct" << std::endl;
	}
	if (outcome == INCORRECT) {
		nIncorrect++;
		std::cout << "Incorrect" << std::endl;
	}
}