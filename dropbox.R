library(rdrop2)

# regenerate the token when needed
rdrop2::drop_auth(new_user = T)
token <- rdrop2::drop_auth()
saveRDS(token, "curated/token.rds")
saveRDS(token, "tests/testthat/curated/token.rds")
