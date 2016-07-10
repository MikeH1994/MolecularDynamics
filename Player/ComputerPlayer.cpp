#include "ComputerPlayer.h"

ComputerPlayer::ComputerPlayer(float money, std::string name, PlayerType character) : GenericPlayer(money, name) {
	this->character = character;
}

GenericPlayer::Moves ComputerPlayer::getMove(bool print) {
	Moves move = ERR;
	switch (character) {
	case NPC_ALEX:
		move = moveAlex();
		break;
	case NPC_SIDNEY:
		move = moveSid();
		break;
	}
	if (print) {
		printMove(move);
	}
	return move;
}

GenericPlayer::Moves ComputerPlayer::moveAlex() {
	if (getScore() < 21) {
		return HIT;
	}
	else {
		return STICK;
	}
}

GenericPlayer::Moves ComputerPlayer::moveSid() {
	if (getScore() < 11) {
		return HIT;
	}
	else {
		return STICK;
	}
}