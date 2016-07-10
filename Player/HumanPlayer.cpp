#include "HumanPlayer.h"

HumanPlayer::HumanPlayer(float money, std::string name) :GenericPlayer(money, name) {
	this->name = name;
}

GenericPlayer::Moves HumanPlayer::getMove(bool print) {
	Moves move;
	if (getScore() <= 21) {
		std::vector<Moves> allowedMoves = getAvailableMoves();
		std::string str = "Enter move for " + getName();
		str += "\n" + getAvailableMovesStr(allowedMoves);
		move = Helper::getUserMove(allowedMoves, str);
	}
	else {
		move = BUST;
	}
	if (print) {
		printMove(move);
	}
	return move;
}

std::vector<GenericPlayer::Moves> HumanPlayer::getAvailableMoves() {
	std::vector<GenericPlayer::Moves> moves = { QUIT,STICK };
	if (getScore() < 21) {
		moves.push_back(HIT);
		moves.push_back(FOLD);
		if (getMoney() >= getBet()) {
			moves.push_back(DOUBLE);
		}
	}
	return moves;
}

std::string HumanPlayer::getAvailableMovesStr(std::vector<Moves> &allowedMoves) {
	std::string str = "Available moves: ";
	for (unsigned int i = 0; i < allowedMoves.size(); i++) {
		str += " '" + Helper::getMoveStr(allowedMoves[i]) + "'; ";
	}
	return str;
}