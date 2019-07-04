function converted_seq = matching_states_convert(seq, matching) 

    % (Not used) it converts all inferred states in a states path (including redundant ones)
    % into their corresponding true ones.
    %
    % INPUTS
    % ------
    % seq             inferred path
    % matching        Boolean matrix obtained with get_states_correspondence.m
    %
    %
    % OUTPUT
    % ------
    % converted_seq   inferred path with the labels changed.
    %
    %
    % Author: Luis Torada Aguilella
    
    
    for s = 1:5
        ms = find(matching(s,:));
        for s2 = ms
            seq(seq == s2) = s*100;
        end
    end
    converted_seq = seq / 100;
    
end