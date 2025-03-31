use proc_macro::TokenStream;
use proc_macro::Span;
use quote::quote;
use syn::{Ident, parse_macro_input, DeriveInput};
use proc_macro_crate::{crate_name, FoundCrate};

/// Placeholder derive macro templace for development. This does nothing in itself.
#[proc_macro_derive(Placeholder)]
pub fn derive_coord_placeholder(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let _name = &input.ident;
    quote! {

    }.into()
}

/// Default implementation of the GenomeCoordinate Trait.
/// 
/// The struct must contain a field named 'coordinate' and contain a 'Coordinate' struct.
#[proc_macro_derive(Coord)]
pub fn derive_genomic_coord(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;
    let crate_prefix = match crate_name("genome") {
        Ok(FoundCrate::Itself)     =>  Ident::new("crate", Span::call_site().into()),
        Ok(FoundCrate::Name(name)) =>  Ident::new(&name, Span::call_site().into()),
        Err(e)                     =>  panic!("# While parsing CoordEq derive macro: [{e:?}]")
    };

    quote! {
        use #crate_prefix::coordinate::GenomeCoordinate;
        
        impl GenomeCoordinate for #name {
            fn coordinate(&self) -> &'_ Coordinate {
                &self.coordinate
            }
        }
    }.into()
}

/// Default implementation of PartialEq for anything that contains a Coordinate struct.
///
/// Structs are considered equal if their chromosome AND position matches.
#[proc_macro_derive(CoordEq)]
pub fn derive_coord_eq(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;
    quote! {
        // ---- PartialEq Rhs = Self
        impl ::core::cmp::PartialEq<Self> for #name {
            fn eq(&self, other: &Self) -> bool { 
                let (inner, other) = (self.coordinate(), other.coordinate());
                inner.chromosome == other.chromosome && inner.position == other.position
            }
        }
        // ---- Eq
        impl ::core::cmp::Eq for #name {}
    }.into()
}

/// Default implementation of PartialEq between a coordinate struct and a JackknifeBlock
/// 
/// Coordinates are considered equal to the Block, if chromosome matches AND they are contained within the block's position Range.
#[proc_macro_derive(CoordBlockEq)]
pub fn derive_coord_block_eq(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;
    quote! {
        impl PartialEq<crate::jackknife::JackknifeBlock> for #name {
            fn eq(&self, block: &crate::jackknife::JackknifeBlock) -> bool {
                let inner = self.coordinate();
                block.chromosome == inner.chromosome && inner.position >= block.range.start && inner.position < block.range.end
            }
        }
    }.into()
}

#[proc_macro_derive(CoordOrd)]
pub fn derive_coord_ord(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;
    quote! {
        // ---- PartialOrd
        impl ::core::cmp::PartialOrd<Self> for #name {
            fn partial_cmp(&self, other: &Self) -> Option<::core::cmp::Ordering> {
                Some((self.coordinate()).cmp(other.coordinate()))
            }
        }
        // ---- Ord
        impl ::core::cmp::Ord for #name {
            fn cmp(&self, other: &Self) -> ::core::cmp::Ordering {
                let (inner, other) = (self.coordinate(), other.coordinate());
                (inner.chromosome, inner.position).cmp(&(other.chromosome, other.position))
            }
        }
    }.into()
}
/// Default implementation of Hash for any struct implementing the `GenomeCoordinate` trait.
/// 
/// Coordinates are hashed using chromosome and position information.
#[proc_macro_derive(CoordHash)]
pub fn derive_coord_hash(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;
    quote! {
        // ---- Hash
        impl ::core::hash::Hash for #name {
            fn hash<H: ::core::hash::Hasher>(&self, state: &mut H) {
                let inner = self.coordinate();
                inner.chromosome.hash(state);
                inner.position.hash(state);
            }
        }
    }.into()
}

#[proc_macro_derive(CoordBorrow)]
pub fn derive_coord_borrow(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;
    quote! {
        impl ::core::borrow::Borrow<Coordinate> for #name  {
            fn borrow(&self) -> &'_ Coordinate {
                &self.coordinate
            }
        }
    }.into()
}