function split(str, sep) -- splits string by delimiter
	local sep, fields = sep or ":", {}
	local pattern = string.format("([^%s]+)", sep)
	str:gsub(pattern, function(c) fields[#fields+1] = c end)
	return fields
end

function split_transcripts(str) -- splits multiple transcripts entries for a single variant
	local transcripts = {}
	for transcript in str:gmatch("(.-Warnings:)") do
		table.insert(transcripts, transcript)
	end
	return transcripts
end

function abs(val)
	return math.abs(val)
end

function pangolin(entry) -- processes Pangolin scores to check for maximum transcript score
	local transcripts = split_transcripts(entry)
	local maximums = {}
	local max_index = 1
	local max_score = 0
	for i = 1, #transcripts do
		local scores = split(transcripts[i], "|")
		local increase = scores[2]
		local decrease = scores[3]
		local current_score = math.max(math.abs(tonumber(split(increase, ":")[2])), math.abs(tonumber(split(decrease, ":")[2])))
		if current_score > max_score then
			max_index = i
			max_score = current_score
		end
	end
	return transcripts[max_index]
end
