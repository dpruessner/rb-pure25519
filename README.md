# Ruby Curve25519

DJB's [Curve25519](http://cr.yp.to/ecdh.html)  implemented in pure Ruby using
core Bignum for calculations on prime Finite Fields.

## Usage

### Ephermeral DH Exchange

```ruby
require 'rb-pure25519'


skey1 = Rb25519.random_secret_str        # Local
pkey1 = Rb25519.public_key_str(skey1)    # Local, pkey1 shared publically

skey2 = Rb25519.random_secret_str        # Remote
pkey2 = Rb25519.public_key_str(skey2)    # Remote, pkey2 shared publically

shared_secret1 = Rb25519.shared_secret_str(pkey1, skey2)  # Local
shared_secret2 = Rb25519.shared_secret_str(pkey2, skey1)  # Remote

shared_secret1 == shared_secret2   #==> true              # shared_secret used for symmetric crypto

```

## Other notes

There's a number of functions for other Montgomery curves — there's a toy curve in there that I used to 
verify some assumptions in the Montgomery EC ladder and in the number theory
usin Sage Math.

We started using EC25519 when TLS was too heavy for an embedded IoT product: 
[Sprinkl Conserve](http://sprinkl.io).   We were looking for a pure ruby
EC25519.  We landed on using RbNaCl for the server-side (it was available in the
Ubuntu distribution on our cloud platform).  We later returned to the subject to
create this library for anyone who needs a pure ruby EC25519 implementation.




Powered by Curve25519.
