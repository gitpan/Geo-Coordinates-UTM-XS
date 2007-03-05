#!/usr/bin/perl

use Test::More tests => 50001;
BEGIN { use_ok('Geo::Coordinates::UTM::XS') };

use constant maxerror => 1e-3;
sub feq ($$) { abs($_[0] - $_[1]) < maxerror }

open TD, '<', 't/test.dat'
    or die "unable to open test data file 't/test.dat'";

while(<TD>) {
    chomp;
    next if /^\s*(?:#.*)?$/;
    my ($ellipsoid, $latitude, $longitude, $zone, $easting, $northing) = split /\|/;
    my ($z, $e, $n) = latlon_to_utm($ellipsoid, $latitude, $longitude);
    ok($z eq $zone, "zone $.");
    ok(feq($e, $easting), "easting $.");
    ok(feq($n, $northing), "northing $.");
    my ($lon, $lat) = utm_to_latlon($ellipsoid, $zone, $easting, $northing);
    ok(feq($lon, $longitude), "longitude $.");
    ok(feq($lat, $latitude), "latitude $."); 
}

