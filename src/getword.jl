#(get a word)
function getword(buffer, word_count, word_window)
   word_count = word_count + 1
   word = buffer[((length(buffer)-29)-(word_window-1)*30):(length(buffer)-(word_window-1)*30)]
   word_window = word_window + 1
   return word,word_count,word_window
end
