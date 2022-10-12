library(rdrop2)

# regenerate the token when needed
rdrop2::drop_auth(new_user = T)
token <- rdrop2::drop_auth()
saveRDS(token, "token.rds")
saveRDS(token, "tests/testthat/token.rds")
