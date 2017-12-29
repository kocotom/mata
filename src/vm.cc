// TODO: add header

#include <vata2/vm.hh>
#include <vata2/vm-dispatch.hh>


void Vata2::VM::VirtualMachine::run(const Vata2::Parser::Parsed& parsed)
{ // {{{
	for (const auto& parsec : parsed)
	{
		this->run(parsec);
	}
} // run(Parsed) }}}


/**
 * @brief  Executes a parsed section
 *
 * This method is the workhorse of VATA.  It interprets an element of
 * ParsedSection and accordingly modifies the state of the virtual machine.
 *
 * If @p parsec is not of the type CODE, an object of the corresponding type is
 * constructed and saved into the virtual machine.
 *
 * If @p parsec is of the type CODE, the program inside is executed, taking into
 * account the state of the virtual machine, and modifying it.
 *
 * @param[in]  parsec  The parsed section to execute
 */
void Vata2::VM::VirtualMachine::run(const Vata2::Parser::ParsedSection& parsec)
{ // {{{
	if (parsec.type == "CODE") {
		this->run_code(parsec);
		return;
	}

	DEBUG_PRINT_LN(parsec.type << "\n");

	VMDispatcherFunc dispatch = Vata2::VM::find_dispatcher(parsec.type);
	std::string name;
	if (!util::haskey(parsec.dict, "Name")) {
		WARN_PRINT("constructing an entity without a name; "
			"the result will be discarded");
	}
	else {
		const auto& vec_of_strs = parsec.dict.at("Name");
		if (vec_of_strs.size() == 0) {
			throw std::runtime_error("the \"name\" attribute provided without value");
		}
		if (vec_of_strs.size() > 1) {
			throw std::runtime_error(
				"the \"name\" attribute provided with multiple values");
		}
		else {
			name = vec_of_strs[0];
		}
	}

	VMFuncArgs args = {{"Parsec", &parsec}};
	VMValue val = dispatch("construct", args);
	if (!name.empty()) {
		this->mem[name] = val;
	}

	assert(false);
} // run(ParsedSection) }}}


void Vata2::VM::VirtualMachine::run_code(
	const Vata2::Parser::ParsedSection& parsec)
{ // {{{
	DEBUG_PRINT("VATA-CODE START");
	for (const auto& line : parsec.body)
	{
		DEBUG_PRINT(std::to_string(line));
		this->execute_line(line);
	}
	DEBUG_PRINT("VATA-CODE END");
	DEBUG_PRINT("Stack: " + std::to_string(this->exec_stack));
} // run_code(ParsedSection) }}}


void Vata2::VM::VirtualMachine::execute_line(
	const Vata2::Parser::BodyLine& line)
{ // {{{
	for (const std::string& tok : line) {
		this->process_token(tok);
	}
} // execute_line(BodyLine) }}}


void Vata2::VM::VirtualMachine::process_token(
	const std::string& tok)
{ // {{{
	if (")" != tok) { // nothing special
		std::string* str = new std::string(tok);
		VMValue val("string", str);
		this->exec_stack.push(val);
	} else { // closing parenthesis - execute action
		assert(")" == tok);

		std::vector<VMValue> exec_vec;
		bool closed = false;
		while (!this->exec_stack.empty()) {
			const VMValue& st_top = this->exec_stack.top();
			if ("string" == st_top.type) {
				assert(nullptr != st_top.ptr);
				const std::string& val = *static_cast<const std::string*>(st_top.ptr);
				if ("(" == val) {
					closed = true;
					this->exec_stack.pop();
					break;
				}
			}

			exec_vec.insert(exec_vec.begin(), st_top);
			this->exec_stack.pop();
		}

		if (!closed) {
			assert(false);
		}

		this->exec_cmd(exec_vec);
	}
} // process_token(string) }}}


void Vata2::VM::VirtualMachine::exec_cmd(
	const std::vector<VMValue>& exec_vec)
{ // {{{
	DEBUG_PRINT("Executing " + std::to_string(exec_vec));

	// getting the function name
	assert(exec_vec.size() >= 2);
	const VMValue& fnc_val = exec_vec[0];
	assert("string" == fnc_val.type);
	assert(nullptr != fnc_val.ptr);
	const std::string& fnc_name = *static_cast<const std::string*>(fnc_val.ptr);

	// getting the object type (type of the first argument of the function)
	const VMValue& arg1_val = exec_vec[1];
	const std::string& arg1_type = arg1_val.type;

	// constructing the arguments
	std::vector<VMValue> args(exec_vec.begin() + 1, exec_vec.end());

	VMValue ret_val = find_dispatcher(arg1_type)(fnc_name, args);



	// mess
	// VMValue ret_val("result", nullptr);

	this->exec_stack.push(ret_val);
} // exec_cmd(std::vector) }}}
