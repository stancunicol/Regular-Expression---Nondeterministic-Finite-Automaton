#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <stack>
#include <map>
#include <queue>
#include <fstream>

//prioritatea fiecarui operator
int Priority(const char& op) {
    if (op == '*') return 3;
    if (op == '.') return 2;
    if (op == '|') return 1;
    return 0;
}

//functie care verifica daca cuvantul este valid
bool CheckExpression(const std::string& expression) {
    std::stack<uint8_t> characters;
    for (const auto& c : expression) {
        if (Priority(c) == 0)
            characters.push(1);
        else if (Priority(c) == 3) {
            if (characters.empty())
                return false;
        }
        else if (Priority(c) == 2) {
            if (characters.size() < 2)
                return false;
            characters.pop();
        }
        else if (Priority(c) == 1) {
            if (characters.size() < 2)
                return false;
            characters.pop();
        }
    }
    if (characters.size() == 1)
        return true;
    return false;
}

//am creat o functie care adauga concatenarea acolo unde lipseste
std::string AddDotsWhereNeeded(const std::string& expression) {
    std::string result = "";
    for (size_t i = 0; i < expression.size(); i++) {
        char character = expression[i];
        result += character;
        if (i + 1 < expression.size()) {
            char nextCharacter = expression[i + 1];
            if ((isalnum(character) || character == ')' || character == '*') &&
                (isalnum(nextCharacter) || nextCharacter == '('))
                result += '.';
        }
    }
    return result;
}

//creez functia post fixata
std::string PostFixExpression(const std::string& expression) {
    std::string postFixExpression = "";
    std::stack<char> stack;
    for (const char& ch : expression) {
        if (isalnum(ch))
            postFixExpression += ch;
        else if (ch == '(')
            stack.push(ch);
        else if (ch == ')') {
            while (!stack.empty() && stack.top() != '(') {
                postFixExpression += stack.top();
                stack.pop();
            }
            stack.pop();
        }
        else if (ch == '*' || ch == '+' || ch == '.' || ch == '|') {
            while (!stack.empty() && Priority(ch) <= Priority(stack.top())) {
                postFixExpression += stack.top();
                stack.pop();
            }
            stack.push(ch);
        }
    }
    while (!stack.empty()) {
        postFixExpression += stack.top();
        stack.pop();
    }
    return postFixExpression;
}





//clasa pentru automatul finit nedeterminist
class NondeterministicFiniteAutomaton {
protected:
    uint8_t m_start;//starea de start
    uint8_t m_stop;//starea de stop
    std::map<uint8_t, std::map<char, std::set<uint8_t>>> m_transitions;//aici retinem tranzitiile pentru fiecare stare, unde avem stare start -> simbol -> stare stop

public:
    NondeterministicFiniteAutomaton() : m_start{ 0 }, m_stop{ 0 } {}//constructor


    NondeterministicFiniteAutomaton(const NondeterministicFiniteAutomaton& other); //constructor de copiere
    NondeterministicFiniteAutomaton& operator=(const NondeterministicFiniteAutomaton& other);


    void AddTransition(uint8_t from, char symbol, uint8_t to);//aici adaugam tranzitii
    int CreateState();//actualizam contorul
    NondeterministicFiniteAutomaton PostFixToNFA(const std::string& postfix);//creem automatul finit nedeterminist in functie de expresia noastra
    void PrintNFA(const NondeterministicFiniteAutomaton& nfa);//afisam automatul


    std::map<uint8_t, std::map<char, std::set<uint8_t>>> GetTransitions() const {
        return m_transitions;
    }

    uint8_t GetStart() const {
        return m_start;
    }
    uint8_t GetStop() const {
        return m_stop;
    }
};

NondeterministicFiniteAutomaton::NondeterministicFiniteAutomaton(const NondeterministicFiniteAutomaton& other)
    : m_start(other.m_start), m_stop(other.m_stop), m_transitions(other.m_transitions) {
    //std::cout << "Copy constructor called.\n";
}

NondeterministicFiniteAutomaton& NondeterministicFiniteAutomaton::operator=(const NondeterministicFiniteAutomaton& other)
{
    if (this != &other) {
        m_start = other.m_start;
        m_stop = other.m_stop;
        m_transitions = other.m_transitions;
    }
    return *this;
}

void NondeterministicFiniteAutomaton::AddTransition(uint8_t from, char symbol, uint8_t to) {
    m_transitions[from][symbol].insert(to);
}

int NondeterministicFiniteAutomaton::CreateState() {
    static uint8_t stateCount = 0;
    return stateCount++;
}

NondeterministicFiniteAutomaton NondeterministicFiniteAutomaton::PostFixToNFA(const std::string& postfix) {
    std::stack<NondeterministicFiniteAutomaton> nfaStack;

    for (char c : postfix) {
        if (isalpha(c)) {
            NondeterministicFiniteAutomaton nfa;
            nfa.m_start = CreateState();
            nfa.m_stop = CreateState();
            nfa.AddTransition(nfa.m_start, c, nfa.m_stop);
            nfaStack.push(nfa);
        }
        else if (c == '|') {
            NondeterministicFiniteAutomaton nfa2 = nfaStack.top(); nfaStack.pop();
            NondeterministicFiniteAutomaton nfa1 = nfaStack.top(); nfaStack.pop();

            NondeterministicFiniteAutomaton newNFA;
            newNFA.m_start = CreateState();
            newNFA.m_stop = CreateState();

            newNFA.AddTransition(newNFA.m_start, 'L', nfa1.m_start);
            newNFA.AddTransition(newNFA.m_start, 'L', nfa2.m_start);
            newNFA.AddTransition(nfa1.m_stop, 'L', newNFA.m_stop);
            newNFA.AddTransition(nfa2.m_stop, 'L', newNFA.m_stop);

            newNFA.m_transitions.insert(nfa1.m_transitions.begin(), nfa1.m_transitions.end());
            newNFA.m_transitions.insert(nfa2.m_transitions.begin(), nfa2.m_transitions.end());

            nfaStack.push(newNFA);
        }
        else if (c == '.') {
            NondeterministicFiniteAutomaton nfa2 = nfaStack.top(); nfaStack.pop();
            NondeterministicFiniteAutomaton nfa1 = nfaStack.top(); nfaStack.pop();

            nfa1.AddTransition(nfa1.m_stop, 'L', nfa2.m_start);
            nfa1.m_stop = nfa2.m_stop;

            nfa1.m_transitions.insert(nfa2.m_transitions.begin(), nfa2.m_transitions.end());

            nfaStack.push(nfa1);
        }
        else if (c == '*') {
            NondeterministicFiniteAutomaton nfa = nfaStack.top(); nfaStack.pop();

            NondeterministicFiniteAutomaton newNFA;
            newNFA.m_start = CreateState();
            newNFA.m_stop = CreateState();

            newNFA.AddTransition(newNFA.m_start, 'L', nfa.m_start);
            newNFA.AddTransition(nfa.m_stop, 'L', newNFA.m_stop);
            newNFA.AddTransition(newNFA.m_start, 'L', newNFA.m_stop);
            newNFA.AddTransition(nfa.m_stop, 'L', nfa.m_start);

            newNFA.m_transitions.insert(nfa.m_transitions.begin(), nfa.m_transitions.end());

            nfaStack.push(newNFA);

        }
    }

    NondeterministicFiniteAutomaton resultNFA = nfaStack.top();

    resultNFA.m_start = 0;

    return resultNFA;
}

void NondeterministicFiniteAutomaton::PrintNFA(const NondeterministicFiniteAutomaton& nfa) {
    std::cout << "\n\nNon-deterministic Finite Automaton\n";
    std::cout << "Start state: " << (int)nfa.m_start << '\n';
    std::cout << "Stop state: " << (int)nfa.m_stop << '\n';
    std::set<int> states;
    for (const auto& trans : nfa.m_transitions) {
        for (const auto& symbol : trans.second) {
            for (int to : symbol.second) {
                states.insert(trans.first);
                states.insert(to);
            }
        }
    }

    std::cout << "States: ";
    for (const int state : states) {
        std::cout << (int)state << " ";
    }
    std::cout << '\n';
    for (const auto& trans : nfa.m_transitions) {
        for (const auto& symbol : trans.second) {
            for (int to : symbol.second) {
                std::cout << (int)trans.first << " -- " << symbol.first << " --> " << (int)to << '\n';
            }
        }
    }
}


//clasa finite daterminist
class DeterministicFiniteAutomaton : public NondeterministicFiniteAutomaton{
private:
    std::set<uint8_t> m_Q; // Multimea starilor DFA
    std::set<char> m_E;    // Alfabetul
    std::map<uint8_t, std::map<char, uint8_t>> m_Transitions; // Tranzitii DFA
    uint8_t q0;            // Starea initiala DFA
    std::set<uint8_t> F;   // Multimea starilor finale DFA

public:
    DeterministicFiniteAutomaton ConvertFromNFA(const NondeterministicFiniteAutomaton& nfa);
    void PrintAutomaton() const;
    void PrintAutomatonF(std::ostream& out) const;

    bool CheckWord(const std::string& word) const;



    std::map<char, uint8_t>& operator[](uint8_t state) {
        return m_Transitions[state];  // Accesează tranzitiile pentru o stare
    }
};

DeterministicFiniteAutomaton DeterministicFiniteAutomaton::ConvertFromNFA(const NondeterministicFiniteAutomaton& nfa) {
    DeterministicFiniteAutomaton dfa;

    // Obține stările inițiale și alfabetul
    std::set<uint8_t> initialClosure;
    std::map<uint8_t, std::map<char, std::set<uint8_t>>>  nfaTransitions = nfa.GetTransitions();

    // Lambda închidere (epsilon-closure) pentru o stare
    auto epsilonClosure = [&nfaTransitions](const std::set<uint8_t>& states) {
        std::set<uint8_t> closure(states);
        std::stack<uint8_t> stack;
        for (uint8_t state : states) stack.push(state);

        while (!stack.empty()) {
            uint8_t top = stack.top();
            stack.pop();
            if (nfaTransitions[top].count('L')) {
                for (uint8_t next : nfaTransitions[top]['L']) {
                    if (!closure.count(next)) {
                        closure.insert(next);
                        stack.push(next);
                    }
                }
            }
        }
        return closure;  // returneaza setul de stari unde se poate ajunge cu lambda
        };

    // Încărcare alfabet, se adauga in dfa toate simbolurile , dif L
    for (const auto& pair : nfaTransitions) {
        uint8_t state = pair.first;
        const auto& transitions = pair.second;
        for (const auto& transition : transitions) {
            char symbol = transition.first;
            if (symbol != 'L') dfa.m_E.insert(symbol);
        }
    }


    // Creez DFA
    std::map<std::set<uint8_t>, uint8_t> dfaStatesMap;
    std::queue<std::set<uint8_t>> unprocessedStates;

    uint8_t startNFA = nfa.GetStart();
    initialClosure = epsilonClosure({startNFA});  // inchiderea starii initiale
    dfaStatesMap[initialClosure] = dfa.q0 = 0; //starea initiala se not cu 0
    unprocessedStates.push(initialClosure);  //adauga in coada ca sa fie procesat

    uint8_t newStateId = 1;

    while (!unprocessedStates.empty()) {
        std::set<uint8_t> currentSet = unprocessedStates.front();
        unprocessedStates.pop();

        uint8_t currentDfaState = dfaStatesMap[currentSet];
        for (char symbol : dfa.m_E) {
            std::set<uint8_t> moveSet;
            for (uint8_t state : currentSet) {
                if (nfaTransitions[state].count(symbol)) {
                    moveSet.insert(nfaTransitions[state][symbol].begin(), nfaTransitions[state][symbol].end());
                }
            }
            std::set<uint8_t> nextStateSet = epsilonClosure(moveSet);
            if (!nextStateSet.empty()) {
                if (!dfaStatesMap.count(nextStateSet)) {
                    dfaStatesMap[nextStateSet] = newStateId++;
                    unprocessedStates.push(nextStateSet);
                }
                dfa.m_Transitions[currentDfaState][symbol] = dfaStatesMap[nextStateSet];
            }
        }
    }

    for (const auto& pair : dfaStatesMap) {
        const std::set<uint8_t>& subset = pair.first;
        uint8_t dfaState = pair.second;

        for (uint8_t state : subset) {
            uint8_t stopNFA = nfa.GetStop();
            if (state == stopNFA) {
                dfa.F.insert(dfaState);
                break;  // Oprim bucla interioară dacă am găsit starea de stop
            }
        }
    }

    return dfa;
}


void DeterministicFiniteAutomaton::PrintAutomaton() const {
    std::cout << "\n\nDeterministic Finite Automaton\n";
    std::cout << "Start state: " << (int)q0 << '\n';
    std::cout << "Final states: ";
    for (uint8_t f : F) std::cout << (int)f << " ";
    std::cout << "\nTransitions:\n";
    
    for (const auto& pair : m_Transitions) {
        uint8_t state = pair.first;
        const auto& transitions = pair.second;
        for (const auto& transition : transitions) {
            char symbol = transition.first;
            uint8_t nextState = transition.second;
            std::cout << (int)state << " -- " << symbol << " --> " << (int)nextState << '\n';

        }
    }
}

void DeterministicFiniteAutomaton::PrintAutomatonF(std::ostream& out) const {
    out << "\n\nDeterministic Finite Automaton\n";
    out << "Start state: " << (int)q0 << '\n';
    out << "Final states: ";
    for (uint8_t f : F) out << (int)f << " ";
    out << "\nTransitions:\n";

    for (const auto& pair : m_Transitions) {
        uint8_t state = pair.first;
        const auto& transitions = pair.second;
        for (const auto& transition : transitions) {
            char symbol = transition.first;
            uint8_t nextState = transition.second;
            out << (int)state << " -- " << symbol << " --> " << (int)nextState << '\n';
        }
    }

    std::cout << "Automaton printed to file." << std::endl; // Debugging message
}



bool DeterministicFiniteAutomaton::CheckWord(const std::string& word) const {
    uint8_t currentState = q0;
    for (char ch : word) {
        // Accesarea tranzitiei folosind operatorul []
        auto& transitions = m_Transitions.at(currentState);  // Folosește at() pentru a evita adăugarea de tranzitii noi
        auto it = transitions.find(ch);

        if (it == transitions.end()) {
            return false;  // Dacă nu găsește simbolul în tranzitii, returnează false
        }

        currentState = it->second;  // Actualizează starea curentă
    }

    // Verifică dacă starea curentă este una finală
    return F.count(currentState) > 0;
}




int main() {
    std::ifstream file;
    file.open("File.txt");
    //output.open("Output.txt");
    std::string expression;
    if (file.is_open()) {
        file >> expression;
    }
    std::ofstream output("Output.txt");
    if (!output.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }
    std::cout << "The expresion with dots is: ";
    std::string expressionWithDots = AddDotsWhereNeeded(expression);
    std::cout << expressionWithDots << '\n';
    std::string postFixedExpression = PostFixExpression(expressionWithDots);
    std::cout << "The post fixed expression is: " << postFixedExpression << '\n';
    if (CheckExpression(postFixedExpression))
        std::cout << "The expression is valid." << '\n';
    else {
        std::cout << "The expression is not valid." << '\n';
        return 0;
    }


    NondeterministicFiniteAutomaton nfa = nfa.PostFixToNFA(postFixedExpression);
    //nfa.PrintNFA(nfa);
    file.close();



    DeterministicFiniteAutomaton dfa;
    dfa = dfa.ConvertFromNFA(nfa);
    dfa.PrintAutomaton();
    
    dfa.PrintAutomatonF(output);
    output.close();

    std::string word;
    std::cout << "Enter a word to test: ";
    std::cin >> word;
    if (dfa.CheckWord(word)) std::cout << "Word is accepted.\n";
    else std::cout << "Word is rejected.\n";

    return 0;
}


/*
a.b.a.(a.a|b.b)*.c.(a.b)*
a(b|c)*d
*/