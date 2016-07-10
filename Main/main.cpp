#include <vector>
#include "Deck.h"
#include "HiLo.h"
#include "Menu.h"
#include "HumanPlayer.h"

const std::vector<int>  Card::CARDS_VAL_ARRAY = { 11,2,3,4,5,6,7,8,9,10,10,10,10 };
const std::vector<std::string> Card::CARDS_STR_ARRAY = { "Ace","Two", "Three","Four","Five", "Six","Seven","Eight","Nine","Ten","Jack","Queen","King" };
const std::vector<std::string> Card::SUIT_STR_ARRAY = { "Hearts","Diamonds","Clubs","Spades" };

int main() {
	Menu* menu = new Menu();
}