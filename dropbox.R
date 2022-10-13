library(rdrop2)

# regenerate the token when needed
# follow this issue for creating long lived tokens: 
# https://github.com/karthik/rdrop2/issues/201
rdrop2::drop_auth(new_user = T)
token <- rdrop2::drop_auth()
saveRDS(token, file.path("inst", "token.rds"))
