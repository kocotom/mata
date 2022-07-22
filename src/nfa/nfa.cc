/* nfa.cc -- operations for NFA
 *
 * Copyright (c) 2018 Ondrej Lengal <ondra.lengal@gmail.com>
 *
 * This file is a part of libmata.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <algorithm>
#include <list>
#include <unordered_set>

// MATA headers
#include <mata/nfa.hh>
#include <simlib/explicit_lts.hh>

using std::tie;

using namespace Mata::util;
using namespace Mata::Nfa;
using Mata::Nfa::Symbol;

const std::string Mata::Nfa::TYPE_NFA = "NFA";

namespace {
//searches for every state in the set of final states
    template<class T>
    bool contains_final(const T &states, const Nfa &automaton) {
        auto finalState = find_if(states.begin(), states.end(),
                                  [&automaton](State s) { return automaton.has_final(s); });
        return (finalState != states.end());
    }

//checks whether the set has someone with the isFinal flag on
    template<class T>
    bool contains_final(const T &states, const std::vector<bool> &isFinal) {
        return any_of(states.begin(), states.end(), [&isFinal](State s) { return isFinal[s]; });
    }

    void uniont_to_left(StateSet &receivingSet, const StateSet &addedSet) {
        receivingSet.insert(addedSet);
    }

    Simlib::Util::BinaryRelation compute_fw_direct_simulation(const Nfa& aut) {
        Simlib::ExplicitLTS LTSforSimulation;
        Symbol maxSymbol = 0;
        const size_t state_num = aut.get_num_of_states();

        for (State stateFrom = 0; stateFrom < state_num; ++stateFrom) {
            for (const TransSymbolStates &t : aut.get_transitions_from_state(stateFrom)) {
                for (State stateTo : t.states_to) {
                    LTSforSimulation.add_transition(stateFrom, t.symbol, stateTo);
                }
                if (t.symbol > maxSymbol) {
                    maxSymbol = t.symbol;
                }
            }
        }

        // final states cannot be simulated by nonfinal -> we add new selfloops over final states with new symbol in LTS
        for (State finalState : aut.finalstates) {
            LTSforSimulation.add_transition(finalState, maxSymbol + 1, finalState);
        }

        LTSforSimulation.init();
        return LTSforSimulation.compute_simulation();
    }
}

std::ostream &std::operator<<(std::ostream &os, const Mata::Nfa::Trans &trans) { // {{{
    std::string result = "(" + std::to_string(trans.src) + ", " +
                         std::to_string(trans.symb) + ", " + std::to_string(trans.tgt) + ")";
    return os << result;
} // operator<<(ostream, Trans) }}}

////// Alphabet related functions

Symbol OnTheFlyAlphabet::translate_symb(const std::string& str)
{ // {{{
    auto it_insert_pair = symbol_map->insert({str, cnt_symbol});
    if (it_insert_pair.second) { return cnt_symbol++; }
    else { return it_insert_pair.first->second; }
} // OnTheFlyAlphabet::translate_symb }}}

std::list<Symbol> OnTheFlyAlphabet::get_symbols() const
{ // {{{
	std::list<Symbol> result;
	for (const auto& str_sym_pair : *(this->symbol_map))
	{
		result.push_back(str_sym_pair.second);
	}

	return result;
} // OnTheFlyAlphabet::get_symbols }}}

std::list<Symbol> OnTheFlyAlphabet::get_complement(
        const std::set<Symbol>& syms) const
{ // {{{
    std::list<Symbol> result;

    // TODO: could be optimized
    std::set<Symbol> symbols_alphabet;
    for (const auto& str_sym_pair : *(this->symbol_map))
    {
        symbols_alphabet.insert(str_sym_pair.second);
    }

    std::set_difference(
            symbols_alphabet.begin(), symbols_alphabet.end(),
            syms.begin(), syms.end(),
            std::inserter(result, result.end()));

    return result;
} // OnTheFlyAlphabet::get_complement }}}

std::list<Symbol> CharAlphabet::get_symbols() const
{ // {{{
    std::list<Symbol> result;
    for (size_t i = 0; i < 256; ++i)
    {
        result.push_back(i);
    }

    return result;
} // CharAlphabet::get_symbols }}}

std::list<Symbol> CharAlphabet::get_complement(
        const std::set<Symbol>& syms) const
{ // {{{
    std::list<Symbol> result;

    std::list<Symbol> symb_list = this->get_symbols();

    // TODO: could be optimized
    std::set<Symbol> symbols_alphabet(symb_list.begin(), symb_list.end());

    std::set_difference(
            symbols_alphabet.begin(), symbols_alphabet.end(),
            syms.begin(), syms.end(),
            std::inserter(result, result.end()));

    return result;
} // CharAlphabet::get_complement }}}

std::list<Symbol> EnumAlphabet::get_symbols() const
{ // {{{
    std::list<Symbol> result;
    for (const auto& str_sym_pair : this->symbol_map)
    {
        result.push_back(str_sym_pair.second);
    }

    return result;
} // EnumAlphabet::get_symbols }}}

std::list<Symbol> EnumAlphabet::get_complement(
        const std::set<Symbol>& syms) const
{ // {{{
    std::list<Symbol> result;

    // TODO: could be optimized
    std::set<Symbol> symbols_alphabet;
    for (const auto& str_sym_pair : this->symbol_map)
    {
        symbols_alphabet.insert(str_sym_pair.second);
    }

    std::set_difference(
            symbols_alphabet.begin(), symbols_alphabet.end(),
            syms.begin(), syms.end(),
            std::inserter(result, result.end()));

    return result;
} // EnumAlphabet::get_complement }}}

///// Nfa structure related methods

void Nfa::add_trans(State stateFrom, Symbol symbolOnTransition, State stateTo) {
    // TODO: define own exceptions
    if (!is_state(stateFrom) || !is_state(stateTo)) {
        throw std::out_of_range(std::to_string(stateFrom) + " or " + std::to_string(stateTo) + " is not a state.");
    }

    auto transitionFromStateIter = transitionrelation[stateFrom].begin();
    for (; transitionFromStateIter != transitionrelation[stateFrom].end(); ++transitionFromStateIter) {
        if (transitionFromStateIter->symbol == symbolOnTransition) {
            transitionFromStateIter->states_to.insert(stateTo);
            return;
        } else if (transitionFromStateIter->symbol > symbolOnTransition) {
            break;
        }
    }

    transitionrelation[stateFrom].insert(transitionFromStateIter, TransSymbolStates(symbolOnTransition, stateTo));
}

void Nfa::remove_trans(State src, Symbol symb, State tgt)
{
    if (!has_trans(src, symb, tgt))
    {
        throw std::invalid_argument(
                "Transition [" + std::to_string(src) + ", " + std::to_string(symb) + ", " +
                std::to_string(tgt) + "] does not exist.");
    }

    auto transitionFromStateIter{ transitionrelation[src].begin() };
    while (transitionFromStateIter != transitionrelation[src].end())
    {
        if (transitionFromStateIter->symbol == symb)
        {
            transitionFromStateIter->states_to.remove(tgt);
            if (transitionFromStateIter->states_to.empty())
            {
                transitionrelation[src].remove(*transitionFromStateIter);
            }

            return;
        }

        ++transitionFromStateIter;
    }
}

State Nfa::add_new_state() {
    transitionrelation.emplace_back();
    return transitionrelation.size() - 1;
}

void Nfa::remove_epsilon(const Symbol epsilon)
{
    *this = Mata::Nfa::remove_epsilon(*this, epsilon);
}

/// General methods for NFA

bool Mata::Nfa::are_state_disjoint(const Nfa& lhs, const Nfa& rhs)
{ // {{{
    // fill lhs_states with all states of lhs
    std::unordered_set<State> lhs_states;
    lhs_states.insert(lhs.initialstates.begin(), lhs.initialstates.end());
    lhs_states.insert(lhs.finalstates.begin(), lhs.finalstates.end());

    for (size_t i = 0; i < lhs.transitionrelation.size(); i++) {
        lhs_states.insert(i);
        for (const auto& symStates : lhs.transitionrelation[i])
        {
            lhs_states.insert(symStates.states_to.begin(), symStates.states_to.end());
        }
    }

    // for every state found in rhs, check its presence in lhs_states
    for (const auto& rhs_st : rhs.initialstates) {
        if (haskey(lhs_states, rhs_st)) { return false; }
    }

    for (const auto& rhs_st : rhs.finalstates) {
        if (haskey(lhs_states, rhs_st)) { return false; }
    }

    for (size_t i = 0; i < lhs.transitionrelation.size(); i++) {
        if (haskey(lhs_states, i))
            return false;
        for (const auto& symState : lhs.transitionrelation[i]) {
            for (const auto& rhState : symState.states_to) {
                if (haskey(lhs_states,rhState)) {
                    return false;
                }
            }
        }
    }

    // no common state found
    return true;
} // are_disjoint }}}

void Mata::Nfa::make_complete(
        Nfa*             aut,
        const Alphabet&  alphabet,
        State            sink_state)
{
    assert(nullptr != aut);

    std::list<State> worklist(aut->initialstates.begin(),
                              aut->initialstates.end());
    std::unordered_set<State> processed(aut->initialstates.begin(),
                                        aut->initialstates.end());
    worklist.push_back(sink_state);
    processed.insert(sink_state);
    while (!worklist.empty())
    {
        State state = *worklist.begin();
        worklist.pop_front();

        std::set<Symbol> used_symbols;
        if (!aut->trans_empty())
        {
            for (const auto &symb_stateset: (*aut)[state]) {
                used_symbols.insert(symb_stateset.symbol);

                const StateSet &stateset = symb_stateset.states_to;
                for (const auto &tgt_state: stateset) {
                    bool inserted;
                    tie(std::ignore, inserted) = processed.insert(tgt_state);
                    if (inserted) { worklist.push_back(tgt_state); }
                }
            }
        }

        auto unused_symbols = alphabet.get_complement(used_symbols);
        for (Symbol symb : unused_symbols)
        {
            aut->add_trans(state, symb, sink_state);
        }
    }
}

void Mata::Nfa::remove_epsilon(Nfa* result, const Nfa& aut, Symbol epsilon)
{
    assert(nullptr != result);

    // cannot use multimap, because it can contain multiple occurrences of (a -> a), (a -> a)
    StateMap<StateSet> eps_closure;

    // TODO: grossly inefficient
    // first we compute the epsilon closure
    for (size_t i=0; i < aut.get_num_of_states(); ++i)
    {
        for (const auto& trans: aut[i])
        { // initialize
            auto it_ins_pair = eps_closure.insert({i, {i}});
            if (trans.symbol == epsilon)
            {
                StateSet &closure = it_ins_pair.first->second;
                closure.insert(trans.states_to);
            }
        }
    }

    bool changed = true;
    while (changed) { // compute the fixpoint
        changed = false;
        for (size_t i=0; i < aut.get_num_of_states(); ++i)
        {
            for (auto const &trans: aut[i])
            {
                if (trans.symbol == epsilon)
                {
                    StateSet &src_eps_cl = eps_closure[i];
                    for (State tgt : trans.states_to)
                    {
                        const StateSet &tgt_eps_cl = eps_closure[tgt];

                        for (State st: tgt_eps_cl)
                        {
                            if (src_eps_cl.count(st) == 0) changed = true;
                            src_eps_cl.insert(st);
                        }
                    }
                }
            }
        }
    }

    // now we construct the automaton without epsilon transitions
    result->initialstates.insert(aut.initialstates);
    result->finalstates.insert(aut.finalstates);
    State max_state{};
    for (const auto& state_closure_pair : eps_closure) { // for every state
        State src_state = state_closure_pair.first;
        for (State eps_cl_state : state_closure_pair.second) { // for every state in its eps cl
            if (aut.has_final(eps_cl_state)) result->make_final(src_state);
            for (const auto& symb_set : aut[eps_cl_state]) {
                if (symb_set.symbol == epsilon) continue;

                // TODO: this could be done more efficiently if we had a better add_trans method
                for (State tgt_state : symb_set.states_to) {
                    max_state = std::max(src_state, tgt_state);
                    if (result->get_num_of_states() < max_state)
                    {
                        result->increase_size_for_state(max_state);
                    }
                    result->add_trans(src_state, symb_set.symbol, tgt_state);
                }
            }
        }
    }
}

void Mata::Nfa::revert(Nfa* result, const Nfa& aut)
{
    assert(nullptr != result);

    if (aut.get_num_of_states() > result->get_num_of_states()) { result->increase_size(aut.get_num_of_states()); }

    result->initialstates = aut.finalstates;
    result->finalstates = aut.initialstates;

    for (size_t i = 0; i < aut.get_num_of_states(); ++i)
    {
        for (const auto& symStates : aut[i])
            for (const State tgt : symStates.states_to)
                result->add_trans(tgt, symStates.symbol, i);
    }
}

bool Mata::Nfa::is_deterministic(const Nfa& aut)
{
    if (aut.initialstates.size() != 1) { return false; }

    if (aut.trans_empty()) { return true; }

    for (size_t i = 0; i < aut.get_num_of_states(); ++i)
    {
        for (const auto& symStates : aut[i])
        {
            if (symStates.states_to.size() != 1) { return false; }
        }
    }

    return true;
}
bool Mata::Nfa::is_complete(const Nfa& aut, const Alphabet& alphabet)
{
    std::list<Symbol> symbs_ls = alphabet.get_symbols();
    std::unordered_set<Symbol> symbs(symbs_ls.cbegin(), symbs_ls.cend());

    // TODO: make a general function for traversal over reachable states that can
    // be shared by other functions?
    std::list<State> worklist(aut.initialstates.begin(),
                              aut.initialstates.end());
    std::unordered_set<State> processed(aut.initialstates.begin(),
                                        aut.initialstates.end());

    while (!worklist.empty())
    {
        State state = *worklist.begin();
        worklist.pop_front();

        size_t n = 0;      // counter of symbols
        if (!aut.trans_empty()) {
            for (const auto &symb_stateset: aut[state]) {
                ++n;
                if (!haskey(symbs, symb_stateset.symbol)) {
                    throw std::runtime_error(std::to_string(__func__) +
                                             ": encountered a symbol that is not in the provided alphabet");
                }

                const StateSet &stateset = symb_stateset.states_to;
                for (const auto &tgt_state: stateset) {
                    bool inserted;
                    tie(std::ignore, inserted) = processed.insert(tgt_state);
                    if (inserted) { worklist.push_back(tgt_state); }
                }
            }
        }

        if (symbs.size() != n) { return false; }
    }

    return true;
}

std::pair<Word, bool> Mata::Nfa::get_word_for_path(const Nfa& aut, const Path& path)
{
    if (path.empty())
    {
        return {{}, true};
    }

    Word word;
    State cur = path[0];
    for (size_t i = 1; i < path.size(); ++i)
    {
        State newSt = path[i];
        bool found = false;

        if (!aut.trans_empty())
        {
            for (const auto &symbolMap: aut.transitionrelation[cur]) {
                for (State st: symbolMap.states_to) {
                    if (st == newSt) {
                        word.push_back(symbolMap.symbol);
                        found = true;
                        break;
                    }
                }

                if (found) { break; }
            }
        }

        if (!found)
        {
            return {{}, false};
        }

        cur = newSt;    // update current state
    }

    return {word, true};
}

bool Mata::Nfa::is_in_lang(const Nfa& aut, const Word& word)
{
    StateSet cur = aut.initialstates;

    for (Symbol sym : word)
    {
        cur = aut.post(cur, sym);
        if (cur.empty()) { return false; }
    }

    return !are_disjoint(cur, aut.finalstates);
}

/// Checks whether the prefix of a string is in the language of an automaton
bool Mata::Nfa::is_prfx_in_lang(const Nfa& aut, const Word& word)
{
    StateSet cur = aut.initialstates;

    for (Symbol sym : word)
    {
        if (!are_disjoint(cur, aut.finalstates)) { return true; }
        cur = aut.post(cur, sym);
        if (cur.empty()) { return false; }
    }

    return !are_disjoint(cur, aut.finalstates);
}

WordSet Mata::Nfa::Nfa::get_shortest_words() const
{
    // Map mapping states to a set of the shortest words accepted by the automaton from the mapped state.
    ShortestWordsMap shortest_words_map{ *this };

    // Get the shortest words for all initial states accepted by the whole automaton (not just a part of the automaton).
    return shortest_words_map.get_shortest_words_for_states(this->initialstates);
}

/// serializes Nfa into a ParsedSection
Mata::Parser::ParsedSection Mata::Nfa::serialize(
        const Nfa&                aut,
        const SymbolToStringMap*  symbol_map,
        const StateToStringMap*   state_map)
{assert(false);}

bool Mata::Nfa::is_lang_empty(const Nfa& aut, Path* cex)
{ // {{{
    std::list<State> worklist(
            aut.initialstates.begin(), aut.initialstates.end());
    std::unordered_set<State> processed(
            aut.initialstates.begin(), aut.initialstates.end());

    // 'paths[s] == t' denotes that state 's' was accessed from state 't',
    // 'paths[s] == s' means that 's' is an initial state
    std::map<State, State> paths;
    for (State s : worklist)
    {	// initialize
        paths[s] = s;
    }

    while (!worklist.empty())
    {
        State state = worklist.front();
        worklist.pop_front();

        if (haskey(aut.finalstates, state))
        {
            // TODO process the CEX
            if (nullptr != cex)
            {
                cex->clear();
                cex->push_back(state);
                while (paths[state] != state)
                {
                    state = paths[state];
                    cex->push_back(state);
                }

                std::reverse(cex->begin(), cex->end());
            }

            return false;
        }

        if (aut.trans_empty())
            continue;

        for (const auto& symb_stateset : aut[state])
        {
            const StateSet& stateset = symb_stateset.states_to;
            for (const auto& tgt_state : stateset)
            {
                bool inserted;
                tie(std::ignore, inserted) = processed.insert(tgt_state);
                if (inserted)
                {
                    worklist.push_back(tgt_state);
                    // also set that tgt_state was accessed from state
                    paths[tgt_state] = state;
                }
                else
                {
                    // the invariant
                    assert(haskey(paths, tgt_state));
                }
            }
        }
    }

    return true;
} // is_lang_empty }}}

bool Mata::Nfa::is_lang_empty_cex(const Nfa& aut, Word* cex)
{
    assert(nullptr != cex);

    Path path = { };
    bool result = is_lang_empty(aut, &path);
    if (result) { return true; }
    bool consistent;
    tie(*cex, consistent) = get_word_for_path(aut, path);
    assert(consistent);

    return false;
}

void Mata::Nfa::complement_in_place(Nfa& aut) {
    StateSet newFinalStates;

    for (State q = 0; q < aut.transitionrelation.size(); ++q) {
        if (aut.finalstates.count(q) == 0) {
            newFinalStates.insert(q);
        }
    }

    aut.finalstates = newFinalStates;
}

void Mata::Nfa::minimize(Nfa *res, const Nfa& aut) {
    //compute the minimal deterministic automaton, Brzozovski algorithm
    Nfa inverted;
    invert(&inverted, aut);
    determinize(res, inverted);
}

void Mata::Nfa::invert(Nfa *invertedAutomaton, const Nfa& aut) {
    const size_t states_num = aut.get_num_of_states();
    for (State i = 0; i < states_num; ++i)
        invertedAutomaton->add_new_state();

    for (State sourceState = 0; sourceState < aut.transitionrelation.size(); ++sourceState)
    {
        for (const TransSymbolStates &transition: aut.transitionrelation[sourceState])
        {
            for(State targetState: transition.states_to)
                invertedAutomaton->add_trans(targetState, transition.symbol, sourceState);
        }
    }

    invertedAutomaton->initialstates = aut.finalstates;
    invertedAutomaton->finalstates = aut.initialstates;
}

void Nfa::print_to_DOT(std::ostream &outputStream) const {
    outputStream << "digraph finiteAutomaton {" << std::endl
                 << "node [shape=circle];" << std::endl;

    for (State finalState: finalstates) {
        outputStream << finalState << " [shape=doublecircle];" << std::endl;
    }

    for (State s = 0; s != transitionrelation.size(); ++s) {
        for (const TransSymbolStates &t: transitionrelation[s]) {
            outputStream << s << " -> {";
            for (State sTo: t.states_to) {
                outputStream << sTo << " ";
            }
            outputStream << "} [label=" << t.symbol << "];" << std::endl;
        }
    }

    outputStream << "node [shape=none, label=\"\"];" << std::endl;
    for (State initialState: initialstates) {
        outputStream << "i" << initialState << " -> " << initialState << ";" << std::endl;
    }

    outputStream << "}" << std::endl;
}

Nfa Nfa::read_from_our_format(std::istream &inputStream) {
    Nfa newNFA;

    // maps names of the states from the input to State
    std::unordered_map<std::string,State> nameToState;
    auto getStateFromName = [&nameToState, &newNFA](const std::string& stateName)->State {
        if (nameToState.count(stateName) > 0) { // if we already seen this name
            return nameToState[stateName];
        } else { // if we have not seen the name yet
            State newState = newNFA.add_new_state();
            nameToState[stateName] = newState;
            return newState;
        }
    };

    // the value is first unused Symbol
    Symbol maxSymbol = 0;

    // maps names of the symbols from the input to Symbol
    std::unordered_map<std::string,Symbol> nameToSymbol;
    auto getSymbolFromName = [&nameToSymbol, &maxSymbol](const std::string& symbolName)->State {
        if (nameToSymbol.count(symbolName) > 0) { // if we already seen this name
            return nameToSymbol[symbolName];
        } else { // if we have not seen the name yet
            // maxSymbol is the new Symbol, we increment it at the end
            nameToSymbol[symbolName] = maxSymbol;
            // we increment maxSymbol, but return the old value (== the new Symbol)
            return maxSymbol++;
        }
    };


    // in definition we keep each "identificator : content" from inputStream that are divided by '@'
    std::string definition;
    // before first @ there should be nothing, we ignore it
    std::getline(inputStream, definition, '@');

    // we read the input where @ is the delimiter
    while (std::getline(inputStream, definition, '@')) {
        if (definition[0] == 'c') {
            continue;
        }

        // remove all whitespaces, see https://stackoverflow.com/questions/83439/remove-spaces-from-stdstring-in-c
        definition.erase(std::remove_if(definition.begin(), definition.end(), [](unsigned char x) { return std::isspace(x); }), definition.end());

        // split definition on :
        std::string identificator = definition.substr(0, definition.find(':'));
        std::string content = definition.substr(definition.find(':') + 1);

        // TODO: define own exception for parsing
        if (identificator == "" || content == "") {
            throw std::runtime_error("Some formula in input is not defined as 'identificator : content'");
        }

        if (identificator[0] == 'k') { // kInitialFormula or kFinalFormula
            // we assume that content looks like "sSTATENAME&sSTATENAME&sSTATENAME&..." where each STATENAME defines initial or final state
            std::istringstream contentStream(content);
            std::string stateName;
            while (std::getline(contentStream, stateName, '&')) {
                // TODO if (stateName[0] != 's') { error }
                State stateToAdd = getStateFromName(stateName.substr(1));
                if (identificator == "kInitialFormula") {
                    newNFA.make_initial(stateToAdd);
                } else if (identificator == "kFinalFormula") {
                    newNFA.make_final(stateToAdd);
                } else {
                    // TODO: define own exception for parsing
                    throw std::runtime_error(std::string("Keywords are kInitialFormula and kFinalFormula but there is ") + identificator);
                }
            }
        } else if (identificator[0] == 's') { // transition
            // assumption: in explicit NFA format, we have transitions "sSTATEFROMNAME:aSYMBOLNAME&sSTATETONAME"

            State stateFrom = getStateFromName(identificator.substr(1));

            // split content on '&'
            std::string symbolName = content.substr(0, content.find('&'));
            std::string stateToName = content.substr(content.find('&') + 1);

            if (symbolName[0] != 'a' || stateToName[0] != 's') {
                // TODO: define own exception for parsing
                throw std::runtime_error("Some transition formula in input is not defined as 'sSTATEFROMNAME:aSYMBOLNAME&sSTATETONAME'");
            }

            Symbol symbol = getSymbolFromName(symbolName.substr(1));
            State stateTo = getStateFromName(stateToName.substr(1));
            newNFA.add_trans(stateFrom, symbol, stateTo);
        } else { // we assume that 'f' formulas are not allowed for explicit NFA
            // TODO: define own exception for parsing
            throw std::runtime_error("On the left side of formula we have to start with keyword or state in input");
        }
    }

    return newNFA;
}

TransSequence Nfa::get_trans_as_sequence()
{
    TransSequence trans_sequence{};

    for (State state_from{ 0 }; state_from < transitionrelation.size(); ++state_from)
    {
        for (const auto& transition_from_state: transitionrelation[state_from])
        {
            for (State state_to: transition_from_state.states_to)
            {
                trans_sequence.push_back(Trans{ state_from, transition_from_state.symbol, state_to });
            }
        }
    }

    return trans_sequence;
}

size_t Nfa::get_num_of_trans() const
{
    size_t num_of_transitions{};

    for (const auto& state_transitions: transitionrelation)
    {
        for (const auto& symbol_transitions: state_transitions)
        {
            num_of_transitions += symbol_transitions.states_to.size();
        }
    }

    return num_of_transitions;
}

void Mata::Nfa::uni(Nfa *unionAutomaton, const Nfa &lhs, const Nfa &rhs) {
    *unionAutomaton = rhs;

    StateMap<State> thisStateToUnionState;
    for (State thisState = 0; thisState < lhs.transitionrelation.size(); ++thisState) {
        thisStateToUnionState[thisState] = unionAutomaton->add_new_state();
    }

    for (State thisInitialState : lhs.initialstates) {
        unionAutomaton->make_initial(thisStateToUnionState[thisInitialState]);
    }

    for (State thisFinalState : lhs.finalstates) {
        unionAutomaton->make_final(thisStateToUnionState[thisFinalState]);
    }

    for (State thisState = 0; thisState < lhs.transitionrelation.size(); ++thisState) {
        State unionState = thisStateToUnionState[thisState];
        for (const TransSymbolStates &transitionFromThisState : lhs.transitionrelation[thisState]) {

            TransSymbolStates transitionFromUnionState(transitionFromThisState.symbol, StateSet{});

            for (State stateTo : transitionFromThisState.states_to) {
                transitionFromUnionState.states_to.insert(thisStateToUnionState[stateTo]);
            }

            unionAutomaton->transitionrelation[unionState].push_back(transitionFromUnionState);
        }
    }
}

void Mata::Nfa::intersection(Nfa *res, const Nfa &lhs, const Nfa &rhs, ProductMap*  prod_map) {
    using StatePair = std::pair<State, State>;

    /* If q is the state of this automaton and p of other, then thisAndOtherStateToIntersectState[q][p] = state in intersect
     * representing (q,p). The hash function computes for each pair (q,p) unique number (q + p*(num of states in this automaton))
     * whose std hash is returned.
     *
     * see https://stackoverflow.com/questions/15719084/how-to-use-lambda-function-as-hash-function-in-unordered-map
     */
    auto hashStatePair = [lhs](const StatePair &sp) {
        return std::hash<unsigned long>()(sp.first + sp.second*lhs.transitionrelation.size());
    };
    // TODO probably remove this since we use prod_map as parameter
    //std::unordered_map<StatePair, State, decltype(hashStatePair)> thisAndOtherStateToIntersectState(10, hashStatePair); // TODO default buckets?

    if (prod_map == nullptr) {
        prod_map = new ProductMap();
    }

    std::vector<StatePair> pairsToProcess;

    for (State thisInitialState : lhs.initialstates) {
        for (State otherInitialState : rhs.initialstates) {
            StatePair thisAndOtherInitialStatePair(thisInitialState, otherInitialState);
            State newIntersectState = res->add_new_state();

            (*prod_map)[thisAndOtherInitialStatePair] = newIntersectState;
            pairsToProcess.push_back(thisAndOtherInitialStatePair);

            res->make_initial(newIntersectState);
            if (lhs.has_final(thisInitialState) && rhs.has_final(otherInitialState))
                res->make_initial(newIntersectState);
            if (lhs.has_final(thisInitialState) && rhs.has_final(otherInitialState)) {
                res->make_final(newIntersectState);
            }
        }
    }

    if (lhs.trans_empty() || rhs.trans_empty())
        return;

    while (!pairsToProcess.empty()) {
        StatePair pairToProcess = pairsToProcess.back();
        pairsToProcess.pop_back();

        // TODO rewrite this

        State intersectState = (*prod_map)[pairToProcess];

        auto thisStateTransitionsIter = lhs.transitionrelation[pairToProcess.first].begin();
        auto otherStateTransitionIter = rhs.transitionrelation[pairToProcess.second].begin();
        auto thisStateTransitionsIterEnd = lhs.transitionrelation[pairToProcess.first].end();
        auto otherStateTransitionIterEnd = rhs.transitionrelation[pairToProcess.second].end();
        // find all transitions that have same symbol for first and the second state in the pairToProcess
        while (thisStateTransitionsIter != thisStateTransitionsIterEnd
               && otherStateTransitionIter != otherStateTransitionIterEnd) {
            // first iterator points to transition with smaller symbol, move it until it is either same or further than second iterator
            if (thisStateTransitionsIter->symbol < otherStateTransitionIter->symbol) {
                while (thisStateTransitionsIter != thisStateTransitionsIterEnd
                       && thisStateTransitionsIter->symbol < otherStateTransitionIter->symbol) {
                    ++thisStateTransitionsIter;
                }
                if (thisStateTransitionsIter == thisStateTransitionsIterEnd) {
                    break;
                }
            } else {
                // second iterator points to transition with smaller symbol, move it until it is either same or further than first iterator
                while (otherStateTransitionIter != otherStateTransitionIterEnd
                       && thisStateTransitionsIter->symbol > otherStateTransitionIter->symbol) {
                    ++otherStateTransitionIter;
                }
                if (otherStateTransitionIter == otherStateTransitionIterEnd) {
                    break;
                }
            }

            // check both iterators point to the transitions with same symbol
            if (thisStateTransitionsIter->symbol == otherStateTransitionIter->symbol) {
                // create transition from the pairToProcess to all pairs between states to which first transition goes and states to which second one goes
                TransSymbolStates intersectTransition(thisStateTransitionsIter->symbol);
                for (State thisStateTo : thisStateTransitionsIter->states_to) {
                    for (State otherStateTo : otherStateTransitionIter->states_to) {
                        StatePair intersectStatePairTo(thisStateTo, otherStateTo);
                        State intersectStateTo;
                        if (prod_map->count(intersectStatePairTo) == 0) {
                            intersectStateTo = res->add_new_state();
                            (*prod_map)[intersectStatePairTo] = intersectStateTo;
                            pairsToProcess.push_back(intersectStatePairTo);

                            if (lhs.has_final(thisStateTo) && rhs.has_final(otherStateTo)) {
                                res->make_final(intersectStateTo);
                            }
                        } else {
                            intersectStateTo = (*prod_map)[intersectStatePairTo];
                        }
                        intersectTransition.states_to.insert(intersectStateTo);
                    }
                }
                res->transitionrelation[intersectState].push_back(intersectTransition);

                ++thisStateTransitionsIter;
                ++otherStateTransitionIter;
            }
        }
    }
}

Simlib::Util::BinaryRelation Mata::Nfa::compute_relation(const Nfa& aut, const StringDict& params) {
    if (!haskey(params, "relation")) {
        throw std::runtime_error(std::to_string(__func__) +
                                 " requires setting the \"relation\" key in the \"params\" argument; "
                                 "received: " + std::to_string(params));
    }
    if (!haskey(params, "direction")) {
        throw std::runtime_error(std::to_string(__func__) +
                                 " requires setting the \"direction\" key in the \"params\" argument; "
                                 "received: " + std::to_string(params));
    }

    const std::string& relation = params.at("relation");
    const std::string& direction = params.at("direction");
    if ("simulation" == relation && direction == "forward") {
        return compute_fw_direct_simulation(aut);
    }
    else {
        throw std::runtime_error(std::to_string(__func__) +
                                 " received an unknown value of the \"relation\" key: " + relation);
    }
}

void Mata::Nfa::determinize(
        Nfa*        result,
        const Nfa&  aut,
        SubsetMap*  subset_map)
{
    //assuming all sets states_to are non-empty
    std::vector<std::pair<State, StateSet>> worklist;
    bool deallocate_subset_map = false;
    if (subset_map == nullptr) {
        subset_map = new SubsetMap();
        deallocate_subset_map = true;
    }
    StateSet S0 =  Mata::Util::OrdVector<State>(aut.initialstates.ToVector());
    State S0id = result->add_new_state();
    result->make_initial(S0id);
    std::vector<bool> isFinal(aut.get_num_of_states(), false);//for fast detection of a final state in a set
    for (const auto &q: aut.finalstates) {
        isFinal[q] = true;
    }
    if (contains_final(S0, isFinal)) {
        result->make_final(S0id);
    }
    worklist.emplace_back(std::make_pair(S0id, S0));

    (*subset_map)[Mata::Util::OrdVector<State>(S0)] = S0id;

    if (aut.trans_empty())
        return;

    while (!worklist.empty()) {
        auto Spair = worklist.back();
        worklist.pop_back();
        StateSet S = Spair.second;
        State Sid = Spair.first;
        if (S.empty()) {
            break;//this should not happen assuming all sets states_to are non empty
        }
        Mata::Nfa::Nfa::state_set_post_iterator iterator(S.ToVector(), aut);

        while (iterator.has_next()) {
            auto symbolTargetPair = iterator.next();
            Symbol currentSymbol = symbolTargetPair.first;
            const StateSet &T = symbolTargetPair.second;
            auto existingTitr = subset_map->find(T);
            State Tid;
            if (existingTitr != subset_map->end()) {
                Tid = existingTitr->second;
            } else {
                Tid = result->add_new_state();
                (*subset_map)[Mata::Util::OrdVector<State>(T)] = Tid;
                if (contains_final(T, isFinal)) {
                    result->make_final(Tid);
                }
                worklist.emplace_back(std::make_pair(Tid, T));
            }
            result->transitionrelation[Sid].push_back(TransSymbolStates(currentSymbol, Tid));
        }
    }

    if (deallocate_subset_map) { delete subset_map; }
}

void Mata::Nfa::construct(
        Nfa*                                 aut,
        const Mata::Parser::ParsedSection&  parsec,
        Alphabet*                            alphabet,
        StringToStateMap*                    state_map)
{ // {{{
    assert(nullptr != aut);
    assert(nullptr != alphabet);

    if (parsec.type != Mata::Nfa::TYPE_NFA) {
        throw std::runtime_error(std::string(__FUNCTION__) + ": expecting type \"" +
                                 Mata::Nfa::TYPE_NFA + "\"");
    }

    bool remove_state_map = false;
    if (nullptr == state_map) {
        state_map = new StringToStateMap();
        remove_state_map = true;
    }

    State cnt_state = 0;

    // a lambda for translating state names to identifiers
    auto get_state_name = [state_map, &cnt_state](const std::string& str) {
        auto it_insert_pair = state_map->insert({str, cnt_state});
        if (it_insert_pair.second) { return cnt_state++; }
        else { return it_insert_pair.first->second; }
    };

    // a lambda for cleanup
    auto clean_up = [&]() {
        if (remove_state_map) { delete state_map; }
    };


    auto it = parsec.dict.find("Initial");
    if (parsec.dict.end() != it)
    {
        for (const auto& str : it->second)
        {
            State state = get_state_name(str);
            aut->initialstates.insert(state);
        }
    }


    it = parsec.dict.find("Final");
    if (parsec.dict.end() != it)
    {
        for (const auto& str : it->second)
        {
            State state = get_state_name(str);
            aut->finalstates.insert(state);
        }
    }

    for (const auto& body_line : parsec.body)
    {
        if (body_line.size() != 3)
        {
            // clean up
            clean_up();

            if (body_line.size() == 2)
            {
                throw std::runtime_error("Epsilon transitions not supported: " +
                                         std::to_string(body_line));
            }
            else
            {
                throw std::runtime_error("Invalid transition: " +
                                         std::to_string(body_line));
            }
        }

        State src_state = get_state_name(body_line[0]);
        Symbol symbol = alphabet->translate_symb(body_line[1]);
        State tgt_state = get_state_name(body_line[2]);

        aut->add_trans(src_state, symbol, tgt_state);
    }

    // do the dishes and take out garbage
    clean_up();
} // construct }}}

void Mata::Nfa::construct(
        Nfa*                                 aut,
        const Mata::Parser::ParsedSection&  parsec,
        StringToSymbolMap*                   symbol_map,
        StringToStateMap*                    state_map)
{ // {{{
    assert(nullptr != aut);

    bool remove_symbol_map = false;
    if (nullptr == symbol_map)
    {
        symbol_map = new StringToSymbolMap();
        remove_symbol_map = true;
    }

    auto release_res = [&](){ if (remove_symbol_map) delete symbol_map; };

    OnTheFlyAlphabet alphabet(symbol_map);

    try
    {
        construct(aut, parsec, &alphabet, state_map);
    }
    catch (std::exception&)
    {
        release_res();
        throw;
    }

    release_res();
} // construct(StringToSymbolMap) }}}

bool Nfa::state_set_post_iterator::has_next() const
{
    return (size > 0);
}

Nfa::state_set_post_iterator::state_set_post_iterator(std::vector<State> states, const Nfa &aut): state_vector(std::move(states)), automaton(aut), size(state_vector.size())
{
    transition_iterators.reserve(size);
    min_symbol = limits.maxSymbol;
    for (size_t i=0; i < size;) {
        State q = state_vector[i];
        if (automaton.transitionrelation[q].empty()) { // we keep in state_vector only states that have some transitions
            transition_iterators[i] = transition_iterators[size-1];
            state_vector[i] = state_vector[size-1];
            size--;
            //then process the same i in the next iteration (no i++ here)
        } else {
            transition_iterators[i] = automaton.transitionrelation[q].begin();
            Symbol transitionSymbol =  transition_iterators[i]->symbol;
            if (transitionSymbol < min_symbol) {
                min_symbol = transitionSymbol;
            }
            i++;
        }
    }
}

Nfa::const_iterator Nfa::const_iterator::for_begin(const Nfa* nfa)
{ // {{{
    assert(nullptr != nfa);

    const_iterator result;
    if (nfa->nothing_in_trans())
    {
        result.is_end = true;
        return result;
    }

    result.nfa = nfa;
    result.trIt = 0;
    assert(!nfa->get_transitions_from_state(0).empty());
    result.tlIt = nfa->get_transitions_from_state(0).begin();
    assert(!nfa->get_transitions_from_state(0).begin()->states_to.empty());
    result.ssIt = result.tlIt->states_to.begin();

    result.refresh_trans();

    return result;
} // for_begin }}}

Nfa::const_iterator Nfa::const_iterator::for_end(const Nfa* /* nfa*/)
{ // {{{
    const_iterator result;
    result.is_end = true;
    return result;
} // for_end }}}

Nfa::const_iterator& Nfa::const_iterator::operator++()
{ // {{{
    assert(nullptr != nfa);

    ++(this->ssIt);
    const StateSet& state_set = this->tlIt->states_to;
    assert(!state_set.empty());
    if (this->ssIt != state_set.end())
    {
        this->refresh_trans();
        return *this;
    }

    // out of state set
    ++(this->tlIt);
    const TransitionList& tlist = this->nfa->get_transitions_from_state(this->trIt);
    assert(!tlist.empty());
    if (this->tlIt != tlist.end())
    {
        this->ssIt = this->tlIt->states_to.begin();

        this->refresh_trans();
        return *this;
    }

    // out of transition list
    ++this->trIt;
    assert(!this->nfa->nothing_in_trans());
    while (this->trIt < this->nfa->transitionrelation.size() &&
        this->nfa->get_transitions_from_state(this->trIt).empty())
    {
        ++this->trIt;
    }

    if (this->trIt < this->nfa->transitionrelation.size())
    {
        this->tlIt = this->nfa->get_transitions_from_state(this->trIt).begin();
        assert(!this->nfa->get_transitions_from_state(this->trIt).empty());
        const StateSet& new_state_set = this->tlIt->states_to;
        assert(!new_state_set.empty());
        this->ssIt = new_state_set.begin();

        this->refresh_trans();
        return *this;
    }

    // out of transitions
    this->is_end = true;

    return *this;
} // operator++ }}}

/**
 * Returns the min_symbol and its post (subset construction), advances the min_symbol to the next minimal symbol,
 * If it meets a state with no more transitions, it swaps it with the last state and decreases the size.
 * size == 0 means no more post.
 * It would be nice to make this parameterized by a functor that defines what is to be done with the individual posts
 * Alternatively, one might return a vector of references to those individual posts, looks like a good idea actually, the cost of this vector should be small enough
 */
 std::pair<Symbol,const StateSet> Nfa::state_set_post_iterator::next() {
    StateSet post;
    Symbol newMinSymbol = limits.maxSymbol;
    for (size_t i=0; i < size;) {
        Symbol transitionSymbol = transition_iterators[i]->symbol;
        if (transitionSymbol == min_symbol) {
            uniont_to_left(post, transition_iterators[i]->states_to);
            transition_iterators[i]++;
            if (transition_iterators[i] != automaton.transitionrelation[state_vector[i]].end()) {
                Symbol nextTransitionSymbol = transition_iterators[i]->symbol;
                if (nextTransitionSymbol < newMinSymbol) {
                    newMinSymbol = nextTransitionSymbol;
                }
                i++;
            } else {
                //swap the ith state with the last state and decrease the size
                transition_iterators[i] = transition_iterators[size-1];
                state_vector[i] = state_vector[size-1];
                size--;
                //do the same i again, the previously last state is there now (no i++)
            }
        } else {
            if (transitionSymbol < newMinSymbol) {
                newMinSymbol = transitionSymbol;
            }
            i++;
        }
    }
    std::pair<Symbol,StateSet> result(min_symbol, post);
    min_symbol = newMinSymbol;
    return result;
}

std::ostream& Mata::Nfa::operator<<(std::ostream& os, const Nfa& nfa)
{ // {{{
    return os << std::to_string(serialize(nfa));
} // Nfa::operator<<(ostream) }}}

std::ostream& std::operator<<(std::ostream& os, const Mata::Nfa::NfaWrapper& nfa_wrap)
{ // {{{
	os << "{NFA wrapper|NFA: " << nfa_wrap.nfa << "|alphabet: " << nfa_wrap.alphabet <<
		"|state_dict: " << std::to_string(nfa_wrap.state_dict) << "}";
	return os;
} // operator<<(NfaWrapper) }}}

WordSet ShortestWordsMap::get_shortest_words_for_states(const StateSet& states) const
{
    std::set <Word> result{};

    if (!shortest_words_map.empty())
    {
        WordLength shortest_words_length{-1};

        for (State state: states)
        {
            const auto& current_shortest_words_map{shortest_words_map.find(state)};
            if (current_shortest_words_map == shortest_words_map.end()) {
                continue;
            }

            const auto& state_shortest_words_map{current_shortest_words_map->second};
            if (result.empty() || state_shortest_words_map.first < shortest_words_length) // Find a new set of the shortest words.
            {
                result = state_shortest_words_map.second;
                shortest_words_length = state_shortest_words_map.first;
            }
            else if (state_shortest_words_map.first == shortest_words_length)
            {
                // Append the shortest words from other state of the same length to the already found set of the shortest words.
                result.insert(state_shortest_words_map.second.begin(),
                              state_shortest_words_map.second.end());
            }
        }

    }

    return result;
}

void Mata::Nfa::ShortestWordsMap::insert_initial_lengths()
{
    if (!reversed_automaton.initialstates.empty())
    {
        for (State state: reversed_automaton.initialstates)
        {
            shortest_words_map.insert(std::make_pair(state, std::make_pair(0, WordSet{ Word{} })));
        }

        processed.insert(reversed_automaton.initialstates.begin(), reversed_automaton.initialstates.end());
        lifo_queue.insert(lifo_queue.end(), reversed_automaton.initialstates.begin(),
                          reversed_automaton.initialstates.end());
    }
}

void ShortestWordsMap::compute()
{
    State state{};
    while (!lifo_queue.empty())
    {
        state = lifo_queue.front();
        lifo_queue.pop_front();

        // Compute the shortest words for the current state.
        compute_for_state(state);
    }
}

void ShortestWordsMap::compute_for_state(const State state)
{
    const LengthWordsPair& dst{ map_default_shortest_words(state) };
    WordLength dst_length_plus_one{ dst.first + 1 };
    LengthWordsPair act{};

    for (const TransSymbolStates& transition: reversed_automaton.get_transitions_from_state(state))
    {
        for (State state_to: transition.states_to)
        {
            const LengthWordsPair& orig{ map_default_shortest_words(state_to) };
            act = orig;

            if ((act.first == -1) || (dst_length_plus_one < act.first))
            {
                // Found new shortest words after appending transition symbols.
                act.second.clear();
                update_current_words(act, dst, transition.symbol);
            }
            else if (dst_length_plus_one == act.first)
            {
                // Append transition symbol to increase length of the shortest words.
                update_current_words(act, dst, transition.symbol);
            }

            if (orig.second != act.second)
            {
                shortest_words_map[state_to] = act;
            }

            if (processed.find(state_to) == processed.end())
            {
                processed.insert(state_to);
                lifo_queue.push_back(state_to);
            }
        }
    }
}

void ShortestWordsMap::update_current_words(LengthWordsPair& act, const LengthWordsPair& dst, const Symbol symbol)
{
    for (Word word: dst.second)
    {
        word.insert(word.begin(), symbol);
        act.second.insert(word);
    }
    act.first = dst.first + 1;
}

void SegNfa::Segmentation::process_state_depth_pair(StateDepthPair& state_depth_pair,
                                            std::deque<StateDepthPair>& worklist)
{
    auto outgoing_transitions{ automaton.get_transitions_from_state(state_depth_pair.state) };
    for (const auto& state_transitions: outgoing_transitions)
    {
        if (state_transitions.symbol == epsilon)
        {
            handle_epsilon_transitions(state_depth_pair, state_transitions, worklist);
        }
        else // Handle other transitions.
        {
            add_transitions_to_worklist(state_transitions, state_depth_pair.depth, worklist);
        }
    }
}

void SegNfa::Segmentation::handle_epsilon_transitions(const StateDepthPair& state_depth_pair,
                                              const TransSymbolStates& state_transitions,
                                              std::deque<StateDepthPair>& worklist)
{
    epsilon_depth_transitions.insert(std::make_pair(state_depth_pair.depth, TransSequence{}));
    for (State target_state: state_transitions.states_to)
    {
        epsilon_depth_transitions[state_depth_pair.depth].push_back(
                Trans{ state_depth_pair.state, state_transitions.symbol, target_state }
        );
        worklist.push_back(StateDepthPair{ target_state, state_depth_pair.depth + 1 });
    }
}

void SegNfa::Segmentation::add_transitions_to_worklist(const TransSymbolStates& state_transitions, EpsilonDepth depth,
                                               std::deque<StateDepthPair>& worklist)
{
    for (State target_state: state_transitions.states_to)
    {
        worklist.push_back(StateDepthPair{ target_state, depth });
    }
}

std::deque<SegNfa::Segmentation::StateDepthPair> SegNfa::Segmentation::initialize_worklist() const
{
    std::deque<StateDepthPair> worklist{};
    for (State state: automaton.initialstates)
    {
        worklist.push_back(StateDepthPair{ state, 0 });
    }
    return worklist;
}

StateMap<bool> SegNfa::Segmentation::initialize_visited_map() const
{
    StateMap<bool> visited{};
    const size_t state_num = automaton.get_num_of_states();
    for (State state{ 0 }; state < state_num; ++state)
    {
        visited[state] = false;
    }
    return visited;
}

void SegNfa::Segmentation::split_aut_into_segments()
{
    segments = AutSequence{ epsilon_depth_transitions.size() + 1, automaton };

    // Construct segment automata.
    std::unique_ptr<const TransSequence> depth_transitions{};
    for (size_t depth{ 0 }; depth < epsilon_depth_transitions.size(); ++depth)
    {
        // Split the left segment from automaton into a new segment.
        depth_transitions = std::make_unique<const TransSequence>(epsilon_depth_transitions[depth]);
        for (const auto& transition: *depth_transitions)
        {
            update_current_segment(depth, transition);
            propagate_to_other_segments(depth, transition);
        }
    }

    trim_segments();
}

void SegNfa::Segmentation::trim_segments()
{
    for (auto& seg_aut: segments)
    {
        seg_aut.trim();
        seg_aut.remove_epsilon(epsilon);
    }
}

void SegNfa::Segmentation::update_current_segment(const size_t current_depth, const Trans& transition)
{
    assert(transition.symb == epsilon);
    assert(segments[current_depth].has_trans(transition));

    segments[current_depth].set_final(transition.src);
    segments[current_depth].remove_trans(transition);
}

void SegNfa::Segmentation::propagate_to_other_segments(const size_t current_depth, const Trans& transition)
{
    const size_t segments_size{ segments.size() };
    for (size_t other_segment_depth{ current_depth + 1 };
         other_segment_depth < segments_size;
         ++other_segment_depth)
    {
        segments[other_segment_depth].remove_trans(transition);
        segments[other_segment_depth].set_initial(transition.tgt);
    }
}

const AutSequence& SegNfa::Segmentation::get_segments()
{
    if (segments.empty()) { split_aut_into_segments(); }

    return segments;
}

void SegNfa::Segmentation::compute_epsilon_depths()
{
    StateMap<bool> visited{ initialize_visited_map() };
    std::deque<StateDepthPair> worklist{ initialize_worklist() };

    while (!worklist.empty())
    {
        StateDepthPair state_depth_pair{ worklist.front() };
        worklist.pop_front();

        if (!visited[state_depth_pair.state])
        {
            visited[state_depth_pair.state] = true;
            process_state_depth_pair(state_depth_pair, worklist);
        }
    }
}
