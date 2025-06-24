use super::IndividualId;

use slotmap::KeyData;

// ------------------------------------------------------------------ //
// ----- Relationship
slotmap::new_key_type! {
    pub struct RelationshipId;
}

impl RelationshipId {
    pub fn to_u64(self) -> u64 {
        self.0.as_ffi()
    }

    pub fn from_u64(id: u64) -> Self {
        RelationshipId::from(KeyData::from_ffi(id))
    }
}

#[derive(Debug, Clone)]
pub struct Relationship {
    pub id: RelationshipId,
    pub from: IndividualId,
    pub to: IndividualId,
}

impl Relationship {
    pub fn new(id: RelationshipId, from: IndividualId, to: IndividualId) -> Self {
        Self { id, from, to }
    }
}
