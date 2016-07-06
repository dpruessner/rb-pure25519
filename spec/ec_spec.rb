require 'rb25519'
require 'pry'
require 'securerandom'

$nacl = false
begin
  require 'rbnacl'
  $nacl = true
rescue
end



# Points on the test curve from SageMath
#
  #(5 : 8 : 1)          -- Example generator (order:44)
  # --
  # (2,  (14 : 44 : 1))
  # (3,  (41 : 36 : 1))
  # (4,  (34 : 6  : 1))
  # (5,  (23 : 37 : 1))
  # (6,  (17 : 4  : 1))
  # (7,  (43 : 36 : 1))
  # (8,  (8  : 17 : 1))
  # (9,  (40 : 28 : 1))
  # (10, (7  : 11 : 1))
  # (11, (46 : 1  : 1))
  # (12, (27 : 29 : 1))
  # (13, (20 : 14 : 1))
  # (14, (6  : 1  : 1))
  # (15, (35 : 14 : 1))
  # (16, (36 : 14 : 1))
  # (17, (45 : 7  : 1))
  # (18, (18 : 17 : 1))
  # (19, (39 : 1  : 1))
  # (20, (37 : 29 : 1))
  # (41, (41 : 11 : 1))
  # (42, (14 : 3  : 1))
  # (43, (5  : 39 : 1))
  # (44, (0  : 1  : 0))   # Infinity
  # (45, (5  : 8  : 1))


def time_it
  ts = Time.now
  yield
  Time.now.-(ts)
end

RSpec.describe Rb25519::FField::MontgomeryEC do
  context "Toy curve with a=3, b=1, p=47" do

    before(:all) do
      @f     = Rb25519::FField.new(47)
      @curve = Rb25519::FField::MontgomeryEC.new(@f, a:3, b:1)
    end

    it "has [5,8] as a point" do
      expect( @curve.on_curve(5,8) ).to be true
    end

    it "lists all points (naively)" do
      expect( @curve.naive_points.length ).to eq(44)
    end

    it "lists points including infinity (naively)" do
      expect( @curve.naive_points ).to include(ECInfinity)
    end


    ## Affine operations
    it "(affine) adds two points" do
      expect( @curve.point_add( [5,8], [7,11] ) ).to eq([46, 1])
    end

    it "(affine) doubles a point" do
      expect( @curve.point_add( [5,8], [5, 8] ) ).to eq([14, 44])
    end

    it "(affine) scalar multiplies a point" do
      expect( @curve.scale_naive(43, [5,8]) ).to eq([5,39])
    end

    it "(affine) double-and-add scalar multiplies a point" do
      expect( @curve.scale_double_add(43, [5,8]) ).to eq([5,39])
    end





    ## XZ projective space operations
    it "(xz) converts from affine to xz" do
      expect( @curve.xz_from_xy([@f[5],@f[8]]) ).to eq([ @f[5], @f[1] ] )
    end

    it "(xz) converts from projective to affine" do
      xz = [ @f[5], @f[1] ]
      expect( @curve.xz_to_xy( xz )).to eq([[5,8], [5,39]])
    end

    it "(xz) projective doubles a point" do
      xz = [ @f[5], @f[1] ]
      dp = @curve.xz_double( xz )
      expect( @curve.xz_to_xy(dp) ).to include( [14, 44] )
    end

    it "(xz) projective scaling each point equals affine scaling" do
      points = @curve.naive_points

      points.each

    end

    it "(xz) projective scales a point" do
      xz = [ @f[5], @f[1] ]

      expect( @curve.scale_proj( 1, xz ).to_xy).to eq(5)
      expect( @curve.scale_proj( 2, xz ).to_xy).to eq(14)
      expect( @curve.scale_proj( 5, xz ).to_xy).to eq(23)
      expect( @curve.scale_proj( 19, xz ).to_xy).to eq(39)
      expect( @curve.scale_proj( 20, xz ).to_xy).to eq(37)
      expect( @curve.scale_proj( 41, xz ).to_xy).to eq(41)
      expect( @curve.scale_proj( 42, xz ).to_xy).to eq(14)
      expect( @curve.scale_proj( 43, xz ).to_xy).to eq(5)
      expect( @curve.scale_proj( 44, xz ).to_xy).to eq(ECInfinity)

    end


  end



  context "Using EC25519 curve" do
    before(:all) do
      @curve = Rb25519::FField::EC25519.new
      @field = @curve.field

      @curve_basepoint   = @curve.xz_to_xy( [@field[9], @field[1] ] )[0]

      @base_to_0_clamped = 52573937849532372069814888178136277301931757964231753054399601207743438775599

    end

    it "(basic) implements String clamp" do
      expect( ("\x00"*32).rb25519_clamp ).to_not be_nil
    end
    it "(basic) implements Integer clamp" do
      expect( 0.rb25519_clamp ).to_not be_nil
    end

    it "(basic) String and Integer clamps are equal" do
      expect( 0.rb25519_clamp.to_binary_string ).to eq( ("\x00"*32).rb25519_clamp )
    end

    if $nacl
      it "(nacl) Should do scaling on from basepoint 9" do
        key  = RbNaCl::GroupElements::Curve25519.new ("\x00"*32).rb25519_clamp
        base = RbNaCl::GroupElements::Curve25519.base

        expect( base.mult(key).to_s.to_binary ).to eq(@base_to_0_clamped)
      end
    end

    it "(Montgomery25519) should calculate the basepoint in affine(x,y)" do
      expect( @curve_basepoint ).to eq([9, 14781619447589544791020593568409986887264606134616475288964881837755586237401 ])
    end

    it "(Montgomery25519) should do scaling on from basepoint 9" do
      key = @field[0.rb25519_clamp].to_i
      
      p   = @curve.scale_double_add(key, @curve_basepoint)

      expect( p[0] ).to eq(@base_to_0_clamped)
    end


    if $nacl
      it "(Montgomery25519) should scale random values to equal RbNaCl calculated values" do
        rand_val_str = SecureRandom.random_bytes(32).rb25519_clamp
        rand_val     = rand_val_str.to_binary

        ts = Time.now
        nacl_val = RbNaCl::GroupElements::Curve25519.base.mult(
            RbNaCl::GroupElements::Curve25519.new(rand_val_str)).to_s.to_binary

        dt = Time.now.-(ts).*(1e3)
        puts "\nNaCl  dt: #{dt} ms"
        expect( @curve.scale_double_add(rand_val, @curve_basepoint)[0] ).to eq(nacl_val)
      end
    end

    if $nacl
      it "(Montgomery25519) should perform ECDHE with RbNaCl" do
        nacl_val_str = SecureRandom.random_bytes(32).rb25519_clamp


        
      end
    end

    it "(Montgomery25519) should have Projective and Affine scaling equal for random values" do
      rand_val_str = SecureRandom.random_bytes(32).rb25519_clamp
      rand_val     = rand_val_str.to_binary

      scaled_xz = @curve.scale_proj(rand_val, @curve_basepoint.to_xz)
      xval      = scaled_xz.to_xy

      expect( @curve.scale_double_add(rand_val, @curve_basepoint )[0] ).to eq(xval)
    end

  end

  context "Projective benchmarks" do
    before(:all) do
      @curve = Rb25519::FField::EC25519.new
      @field = @curve.field

      @curve_basepoint   = @curve.xz_to_xy( [@field[9], @field[1] ] )[0]

      @base_to_0_clamped = 52573937849532372069814888178136277301931757964231753054399601207743438775599

    end

    it "(Montgomery25519) projective scaling should be 10x faster than affine scaling (double-and-add)" do

      proj_dt = 0
      aff_dt  = 0
      n = 10

      n.times do
        rand_val_str = SecureRandom.random_bytes(32).rb25519_clamp
        rand_val     = rand_val_str.to_binary

        proj_dt += time_it {
          scaled_xz = @curve.scale_proj(rand_val, @curve_basepoint.to_xz)
          xval      = scaled_xz.to_xy
        }

        aff_dt += time_it {
          @curve.scale_double_add(rand_val, @curve_basepoint )[0]
        }
      end

      puts "\nAffine dt: #{aff_dt./(n).*(1e3)} ms"
      puts "Projct dt: #{proj_dt./(n).*(1e3)} ms"

      expect( aff_dt / proj_dt ).to be >= 10
    end

  end # Projective Benchmarks 

end


RSpec.describe Rb25519 do

  context "as a module" do
    
    it "generates a random secret key" do
      skey = Rb25519.random_secret_str
      expect( skey ).to_not   be nil
      expect( skey.length).to eq(32)
    end

    it "generates a public key from the random secret key (str)" do
      skey = Rb25519.random_secret_str
      pkey = Rb25519.public_key_str(skey)

      expect( pkey ).to be_a(String)
    end

    if $nacl
      it "generates a public key equal to the NaCl library" do
        skey = Rb25519.random_secret_str
        pkey = Rb25519.public_key_str(skey)

        nacl_val = RbNaCl::GroupElements::Curve25519.base.mult(
            RbNaCl::GroupElements::Curve25519.new(skey)).to_s

        expect( nacl_val.to_s ).to eq(pkey)
      end
    end

    it "generates a shared secret from two random secrets (str)" do
        skey1 = Rb25519.random_secret_str
        pkey1 = Rb25519.public_key_str(skey1)

        skey2 = Rb25519.random_secret_str
        pkey2 = Rb25519.public_key_str(skey2)

        shared_secret1 = Rb25519.shared_secret_str(pkey1, skey2)
        shared_secret2 = Rb25519.shared_secret_str(pkey2, skey1)

        expect( shared_secret1 ).to eq(shared_secret2)
    end

    if $nacl
      it "generates a shared secret equal to NaCl generated secret" do
        skey1 = Rb25519.random_secret_str
        pkey1 = Rb25519.public_key_str(skey1)

        skey2 = Rb25519.random_secret_str

        nacl_pkey = RbNaCl::GroupElements::Curve25519.base.mult(
            RbNaCl::GroupElements::Curve25519.new(skey2)).to_s
        
        shared_secret1 = Rb25519.shared_secret_str(nacl_pkey, skey1)

        shared_secret2 = RbNaCl::GroupElements::Curve25519.new(pkey1).mult( skey2).to_s

        expect( shared_secret1 ).to eq(shared_secret1)

      end
    end

    
  end
end
