make_pw_scatter <- function(res, highlight = NULL) {
    tryCatch(
        expr = pairs(t(res[highlight, ])),
        error = function(e) {
            print("Error in make_pw_scatter:")
            print(e)
        }
    )
}
