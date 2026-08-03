# Two questions that no other check asks, and that both cost a red CI push to
# learn the hard way in the sibling packages.

test_that("every object in the namespace has a help topic", {
  # The question that found seven man pages which existed on Linux and not on
  # Windows, because two topics differing only in case are one file there.
  # Asked here, on the machine where the evidence is destroyed.
  skip_if_not_installed("tools")

  db <- tryCatch(tools::Rd_db("basis7"), error = function(e) NULL)
  skip_if(is.null(db) || !length(db), "package not installed")

  aliases <- unlist(lapply(db, function(rd) {
    tags <- vapply(rd, function(x) attr(x, "Rd_tag"), character(1))
    unlist(lapply(rd[tags == "\\alias"], function(a) trimws(paste(unlist(a),
      collapse = ""
    ))))
  }), use.names = FALSE)

  objs <- getNamespaceExports("basis7")
  missing <- setdiff(objs, aliases)
  expect_identical(missing, character(0))
})


test_that("no two man pages differ only in case", {
  # On a case-insensitive filesystem the second overwrites the first, and
  # nothing complains: R reads the topic from \name{} inside the file, not from
  # the file name. Caught here, where both files can still exist.
  path <- system.file("man", package = "basis7")
  if (!nzchar(path)) path <- "../../man"
  skip_if_not(dir.exists(path), "no man directory")

  files <- list.files(path, pattern = "\\.Rd$")
  skip_if(length(files) == 0L, "no Rd files")
  expect_identical(anyDuplicated(tolower(files)), 0L)
})


test_that("every exported topic is in the pkgdown reference index", {
  # pkgdown enforces this in CI and not locally, so a missing line of YAML is a
  # red push minutes later. Asking here costs a second.
  skip_if_not_installed("pkgdown")
  path <- "../.."
  skip_if_not(file.exists(file.path(path, "_pkgdown.yml")), "not the source tree")
  expect_no_error(pkgdown::check_pkgdown(path))
})
