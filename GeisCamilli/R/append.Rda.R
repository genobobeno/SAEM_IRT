append.Rda <-
function(x, file) {
  old.objects <- load(file, new.env())
  save(list = c(old.objects, deparse(substitute(x))), file = file)
}
