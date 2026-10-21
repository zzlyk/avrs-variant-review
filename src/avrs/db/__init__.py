"""Database session management and schema creation."""
from __future__ import annotations

from pathlib import Path

from sqlalchemy import create_engine, event
from sqlalchemy.orm import Session, sessionmaker

from .models import Base


def make_engine(db_path: str | Path):
    db_path = Path(db_path).expanduser()
    db_path.parent.mkdir(parents=True, exist_ok=True)
    engine = create_engine(f"sqlite:///{db_path}", echo=False)
    # enforce FK constraints in SQLite
    @event.listens_for(engine, "connect")
    def _fk_on(dbapi_conn, _):
        cur = dbapi_conn.cursor()
        cur.execute("PRAGMA foreign_keys=ON")
        cur.execute("PRAGMA journal_mode=WAL")
        cur.close()
    return engine


def init_db(engine) -> None:
    Base.metadata.create_all(engine)


def get_session(engine) -> Session:
    return sessionmaker(bind=engine, expire_on_commit=False)()


def open_db(db_path: str | Path):
    """Convenience: returns (engine, session)."""
    engine = make_engine(db_path)
    init_db(engine)
    return engine, get_session(engine)


def audit(session: Session, actor: str, action: str, entity_type: str | None = None,
          entity_id: str | None = None, details: dict | None = None):
    from .models import AuditLog
    import json as _json
    session.add(AuditLog(
        actor=actor, action=action, entity_type=entity_type, entity_id=str(entity_id) if entity_id else None,
        details_json=_json.dumps(details, ensure_ascii=False) if details else None,
    ))
