using FVSjl
const KEY = "test/fixtures/canonical/snt01.key"
run() = FVSjl.run_keyfile(KEY; variant=FVSjl.Southern(), faithful=true)
run()  # compile
using Profile
Profile.clear_malloc_data()
for _ in 1:8; run(); end
