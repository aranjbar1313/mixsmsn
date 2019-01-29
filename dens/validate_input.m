function outputarg = validate_input (inputarg, default_values)

    inputarg = cell2mat(inputarg);
    input_size = prod(size(inputarg));
    default_size = prod(size(default_values));
    outputarg = inputarg;

    for i = input_size + 1 : default_size

        outputarg(i) = default_values(i);
    end
    outputarg = num2cell(outputarg);
end