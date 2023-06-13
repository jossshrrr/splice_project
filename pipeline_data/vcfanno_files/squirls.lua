function split(str, sep) -- splits string by delimiter
	local sep, fields = sep or ":", {}
	local pattern = string.format("([^%s]+)", sep)
	str:gsub(pattern, function(c) fields[#fields+1] = c end)
	return fields
end

function squirls(entry) -- processes SQUIRLS scores to check for maximum transcript score
    local t = type(entry)
    if t == "string" then
        return entry
    elseif t == "table" then
        local maximum = 0
		local maximum_index = 0
        for i=1,#entry do
            local transcripts = split(entry[i], "|")
            for j=2,#transcripts do
                local scores = split(transcripts[j], "=")
                local score = tonumber(scores[2])
                if score > maximum then
                    maximum = score 
					maximum_index = j
                end
            end
        end 
        return maximum
    end
end

