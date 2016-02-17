function dblArray = symToDouble(symArray)
    dblArray = arrayfun(@double,symArray,'UniformOutput',false);
end