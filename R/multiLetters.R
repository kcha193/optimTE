multiLetters <-
function(x, upper = TRUE) {
    if(upper){
      if (max(x) < 26) {
          LETTERS[x]
      } else {
          newLETTERS = sort(sapply(strsplit(levels(interaction(LETTERS, LETTERS)),
              "\\."), function(x) paste(x, collapse = "", sep = "")))
         newLETTERS[x]
      }
    } else {
      if (max(x) < 26) {
          letters[x]
      } else {
          newLETTERS = sort(sapply(strsplit(levels(interaction(letters, letters)),
              "\\."), function(x) paste(x, collapse = "", sep = "")))
          newLETTERS[x]
      }       
    }
}
